"""Frozen official model inference with verified structure identity and order."""
from pathlib import Path
import sys, json, hashlib, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import MACCSkeys,rdFingerprintGenerator
BASE=Path(__file__).resolve().parent;ROOT=BASE.parent;OFF=ROOT/'models/official_model'
LIB=ROOT/'data/library/public_purchasable_100k.csv'
sys.dont_write_bytecode=True;sys.path.insert(0,str(OFF))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def initialize():
 global LOCK,BUNDLE,REF,STANDARD,DESC,G4
 from src.official import load_model
 LOCK,BUNDLE=load_model()
 from runtime.common import standardize,DESC as descriptors
 STANDARD=standardize;DESC=descriptors
 REF=[]
 for arr in BUNDLE['reference_bits']:
  fp=DataStructs.ExplicitBitVect(len(arr))
  for j in np.flatnonzero(arr):fp.SetBit(int(j))
  REF.append(fp)
 G4=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
def score_chunk(chunk):
 rows=[];xs=[];valid=[]
 for row_id,smi in chunk:
  out={'source_row':row_id,'SMILES':smi,'model_version':LOCK['version'],'status':'invalid_SMILES'}
  s=STANDARD(smi)
  if s is not None and s['heavy_atom_count']>0:
   mol=Chem.MolFromSmiles(s['isomeric_smiles']);fp=MACCSkeys.GenMACCSKeys(mol);bits=np.zeros(167,dtype=np.float32);DataStructs.ConvertToNumpyArray(fp,bits)
   vals=np.array([s[k] for k in DESC],dtype=np.float32)
   similarity=float(np.float32(max(DataStructs.BulkTanimotoSimilarity(G4.GetFingerprint(mol),REF))))
   inside=bool(np.all((vals>=BUNDLE['desc_min'])&(vals<=BUNDLE['desc_max'])))
   out.update(status='OK',standardized_smiles=s['isomeric_smiles'],full_inchikey=s['full_inchikey'],connectivity_key=s['connectivity_key'],murcko_scaffold=s['murcko_scaffold'],max_training_similarity=similarity,descriptors_in_training_range=inside,applicability_domain=bool(similarity>=.4 and inside),**{k:s[k] for k in DESC})
   xs.append(np.r_[bits,vals]);valid.append(len(rows))
  rows.append(out)
 if xs:
  scores=BUNDLE['model'].predict_proba(np.asarray(xs,dtype=np.float32))[:,1]
  for i,p in zip(valid,scores):rows[i].update(ranking_score=float(p),fixed_threshold=LOCK['threshold'],fixed_threshold_class=int(p>=LOCK['threshold']))
 return rows
def main(output_dir, workers=4):
 from src.official import verify_lock,predict
 from src.library_access import read_library
 lock=verify_lock();start=time.perf_counter()
 data,input_receipt=read_library();print('INPUT',input_receipt['input_file'],'rows',len(data),'columns',data.columns.tolist(),flush=True)
 smicol=next(c for c in ['isomeric_smiles','SMILES','smiles','canonical_smiles'] if c in data)
 pairs=list(enumerate(data[smicol].fillna('').astype(str).tolist()))
 chunks=[pairs[i:i+1000] for i in range(0,len(pairs),1000)]
 out=[]
 with ProcessPoolExecutor(max_workers=workers,initializer=initialize) as pool:
  for i,rows in enumerate(pool.map(score_chunk,chunks),1):
   out.extend(rows)
   if i%10==0:print('PROGRESS',i*1000,'elapsed_seconds',round(time.perf_counter()-start,1),flush=True)
 result=pd.DataFrame(out).sort_values('source_row').reset_index(drop=True)
 for c in ['library_id','source_name','source_url','transport_url','availability_status','gpr84_label_status','dataset_version']: 
  if c in data:result[c]=data[c].to_numpy()
 assert len(result)==len(data)==100000
 assert result.status.eq('OK').all() and np.isfinite(result.ranking_score).all()
 # Cross-check the optimized frozen-only feature path against the official public interface.
 test=np.unique(np.r_[np.linspace(0,len(data)-1,16,dtype=int),np.argsort(result.ranking_score.to_numpy())[-8:]])
 reference=predict(data.iloc[test][smicol].tolist())
 assert np.array_equal(reference.ranking_score.to_numpy(),result.iloc[test].ranking_score.to_numpy()),'Score differs from official interface'
 assert np.array_equal(reference.max_training_similarity.to_numpy(),result.iloc[test].max_training_similarity.to_numpy()),'Similarity differs'
 assert np.array_equal(reference.applicability_domain.to_numpy(),result.iloc[test].applicability_domain.to_numpy()),'Domain differs'
 # Preserve the pre-existing blind identities and confirmed training identities as exclusion flags.
 hist=pd.read_csv(ROOT/'data/library/exclusion_flags.csv')
 restricted=set(hist.loc[hist.restricted_holdout_excluded,'library_id'])
 train=pd.read_csv(ROOT/'data/training/training_structures.csv')
 result['restricted_holdout_excluded']=result.library_id.isin(restricted)
 result['training_structure_overlap']=result.full_inchikey.isin(set(train.full_inchikey)) | result.connectivity_key.isin(set(train.connectivity_key))
 result['eligible_for_ranking']=~result.restricted_holdout_excluded & ~result.training_structure_overlap
 ranks=result[result.eligible_for_ranking].sort_values(['ranking_score','library_id'],ascending=[False,True],kind='stable').copy()
 ranks['official_library_rank']=np.arange(1,len(ranks)+1)
 result=result.merge(ranks[['library_id','official_library_rank']],on='library_id',how='left',validate='one_to_one').sort_values('source_row')
 destination=Path(output_dir);destination.mkdir(parents=True,exist_ok=True)
 result.to_parquet(destination/'official_100k_scores.parquet',index=False)
 ranks.head(1000).to_csv(destination/'official_score_top1000.csv',index=False,encoding='utf-8-sig',float_format='%.17g')
 hist_ids=set(pd.read_csv(ROOT/'data/evidence/03_docking_top60.csv').library_id)
 top_ids=set(ranks.head(1000).library_id)
 receipt={'status':'COMPLETE','model_version':lock['version'],'model_sha256':sha(OFF/'model.joblib'),'library_input':input_receipt,'rows':len(result),'valid_rows':int(result.status.eq('OK').sum()),'duplicates_full_inchikey':int(result.full_inchikey.duplicated().sum()),'duplicates_connectivity':int(result.connectivity_key.duplicated().sum()),'positive_at_fixed_threshold':int(result.fixed_threshold_class.sum()),'applicability_domain_rows':int(result.applicability_domain.sum()),'restricted_excluded':int(result.restricted_holdout_excluded.sum()),'training_overlap':int(result.training_structure_overlap.sum()),'eligible_ranking_rows':len(ranks),'interface_score_check_n':len(test),'interface_max_abs_score_delta':0.,'interface_similarity_max_abs_delta':0.,'admet_evidence_available_in_top1000':len(top_ids & set(pd.read_csv(ROOT/'data/evidence/02_admet_top300.csv').library_id)),'docking_evidence_available_in_top1000':len(top_ids & hist_ids),'elapsed_seconds':time.perf_counter()-start,'training_performed':False,'old_results_overwritten':False,'wet_labels_added':False}
 verify_lock()
 (destination/'OFFICIAL_SCREEN_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(receipt,ensure_ascii=False),flush=True)
