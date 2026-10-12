from pathlib import Path
import hashlib,json,datetime
import numpy as np
import pandas as pd
from rdkit import Chem,RDLogger
from rdkit.Chem import Crippen,Descriptors,Lipinski,rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')
OUT=Path(__file__).resolve().parents[1]; ROOT=OUT.parent.parent; DATA=ROOT/'gpr84_ml_dataset_v3'
CONF=json.loads((OUT/'config.json').read_text(encoding='utf8'))
DESC=['molecular_weight','clogp','tpsa','hbd','hba','rotatable_bonds','ring_count','fraction_csp3','formal_charge','heavy_atom_count']
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):return [clean(x) for x in v]
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    return v
def dump(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(clean(v),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def write(d,rel):
    p=OUT/rel;p.parent.mkdir(parents=True,exist_ok=True);d.to_csv(p,index=False,encoding='utf-8-sig')
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def union(v):return '|'.join(sorted({t for x in v if pd.notna(x) for t in str(x).split('|') if t and t!='nan'}))
def txt(x):return '' if pd.isna(x) else str(x).strip()
def num(x):
    try:return float(x)
    except (ValueError,TypeError):return np.nan
def truth(x):return str(x).lower() in ['true','1']
def standardize(s,racemate=False):
    m=Chem.MolFromSmiles(str(s))
    if m is None:return None
    m=rdMolStandardize.Uncharger().uncharge(rdMolStandardize.FragmentParent(m));Chem.SanitizeMol(m)
    if racemate:Chem.RemoveStereochemistry(m)
    key=Chem.MolToInchiKey(m)
    return {'canonical_smiles':Chem.MolToSmiles(m,isomericSmiles=False),'isomeric_smiles':Chem.MolToSmiles(m,isomericSmiles=True),'full_inchikey':key,'connectivity_key':key.split('-')[0],'murcko_scaffold':Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m),isomericSmiles=False),**dict(zip(DESC,[Descriptors.MolWt(m),Crippen.MolLogP(m),rdMolDescriptors.CalcTPSA(m),Lipinski.NumHDonors(m),Lipinski.NumHAcceptors(m),Lipinski.NumRotatableBonds(m),Lipinski.RingCount(m),rdMolDescriptors.CalcFractionCSP3(m),Chem.GetFormalCharge(m),m.GetNumHeavyAtoms()]))}
def evidence_label(r,rawmap):
    reasons=[]
    if r.structure_status!='standardized' or not truth(r.small_molecule_eligible):reasons.append('invalid_structure')
    if txt(r.target_organism).lower() not in ['human','homo sapiens']:reasons.append('non_human')
    if not truth(r.wildtype) or txt(r.variant_mutation):reasons.append('mutant')
    if not np.isfinite(num(r.target_confidence)) or num(r.target_confidence)<8:reasons.append('target_confidence')
    if txt(r.data_validity_comment):reasons.append('validity_flag')
    if truth(r.potential_duplicate):reasons.append('source_duplicate')
    if r.pharmacological_role not in ['antagonist','inverse_agonist_with_measured_inhibition']:reasons.append('not_confirmed_inhibitory_role')
    if reasons:return np.nan,'|'.join(reasons)
    raw=rawmap.get(txt(r.get('chembl_activity_id')).split('.')[0],{})
    desc=txt(r.assay_description).lower()
    direct=txt(raw.get('standard_text_value')).lower() in ['inactive','not active'] and any(x in desc for x in ['receptor activation','beta-arrestin','gtp','camp','calcium mobilization']) and 'antagonist' in desc and truth(r.single_concentration) and num(r.tested_concentration_nm)==10000 and r.manual_review_status=='reviewed_direct_inactive_at_10uM'
    if direct:return 0,'explicit_not_active_at_10uM'
    if txt(r.quality_flag)=='chembl_0_01nm_sentinel':return np.nan,'sentinel'
    if txt(r.endpoint).upper() not in ['IC50','KI','KD','KB']:return np.nan,'unsupported_endpoint'
    value=num(r.standard_value_nm);rel=txt(r.relation) or '='
    if not np.isfinite(value) or value<=0 or r.standard_units!='nM':return np.nan,'invalid_value'
    if rel in ['=','~']:
        if value<=1000:return 1,'exact_active'
        if value>=10000:return 0,'exact_inactive'
        return np.nan,'gray_zone'
    if rel in ['<','<='] and value<=1000:return 1,'bound_active'
    if rel in ['>','>='] and value>=10000:return 0,'bound_inactive'
    return np.nan,'indecisive_bound'
def aggregate(m,structures):
    rows=[];excluded=[]
    for key,g in m[m.record_label.notna() & m.std_key.notna()].groupby('std_key',sort=True):
        conflict=g.record_label.nunique()>1
        exact=g[g.relation.isin(['=','~']) & g.label_basis_record.ne('explicit_not_active_at_10uM')]
        for _,h in exact.groupby(['endpoint','assay_family'],dropna=False):
            p=9-np.log10(pd.to_numeric(h.standard_value_nm,errors='coerce'))
            conflict|=len(p)>1 and p.max()-p.min()>1+1e-10
        if conflict:
            excluded.extend({'measurement_id':r.measurement_id,'reason':'conflicting_binary_or_gt_1log_same_endpoint_family'} for _,r in g.iterrows());continue
        rec=structures[key].copy();rec.update({'compound_id':'V3_'+key,'label':int(g.record_label.iloc[0]),'publication_group':union(g.publication_group_id),'measurement_ids':union(g.measurement_id),'first_publication_year':g.publication_year.min(),'last_publication_year':g.publication_year.max(),'evidence_count':len(g),'label_basis':union(g.label_basis_record),'source_list':union(g.source_list),'merged_compound_ids':union(g.compound_id)})
        rows.append(rec)
    return pd.DataFrame(rows),pd.DataFrame(excluded)
def connected_groups(d,m,scaffold=False):
    keys=list(d.full_inchikey);parent={k:k for k in keys}
    def find(k):
        while parent[k]!=k:parent[k]=parent[parent[k]];k=parent[k]
        return k
    def merge(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    owners={}
    for _,r in d.iterrows():
        tokens=['C:'+r.connectivity_key]
        if scaffold:tokens+=['S:'+txt(r.murcko_scaffold)]
        else:
            pubs=union(m.loc[m.connectivity_key.eq(r.connectivity_key),'publication_group_id']).split('|')
            if not any(pubs):raise ValueError('Unknown publication '+r.compound_id)
            tokens+=['P:'+p for p in pubs if p]
        for t in tokens:
            if t in owners:merge(r.full_inchikey,owners[t])
            else:owners[t]=r.full_inchikey
    return np.array([find(k) for k in keys])
