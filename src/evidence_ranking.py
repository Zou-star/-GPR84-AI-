"""Join frozen molecule-specific downstream evidence to new official scores."""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import FilterCatalog,QED,rdFingerprintGenerator
BASE=Path(__file__).resolve().parent;ROOT=BASE.parent;OFF=ROOT/'models/official_model'
H=ROOT/'data/evidence'
sys.dont_write_bytecode=True;sys.path.insert(0,str(OFF))

def read(p):return pd.read_csv(p,float_precision='round_trip')
def write(d,n):d.to_csv(S/n,index=False,encoding='utf-8-sig',float_format='%.17g')
def main(output_dir):
 global S
 S=Path(output_dir)
 from src.official import predict,verify_lock
 lock=verify_lock();scores=pd.read_parquet(S/'official_100k_scores.parquet')
 ordered=scores[scores.eligible_for_ranking].sort_values(['ranking_score','library_id'],ascending=[False,True]).copy()
 ordered['lipinski_violations']=(ordered.molecular_weight.gt(500).astype(int)+ordered.clogp.gt(5).astype(int)+ordered.hbd.gt(5).astype(int)+ordered.hba.gt(10).astype(int))
 ordered['basic_property_pass']=ordered.molecular_weight.between(150,550)&ordered.clogp.between(-1,5)&ordered.tpsa.le(120)&ordered.rotatable_bonds.le(10)&ordered.lipinski_violations.le(1)
 params=FilterCatalog.FilterCatalogParams();params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS);params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK);catalog=FilterCatalog.FilterCatalog(params)
 chosen=[];evaluated=[]
 for _,r in ordered[ordered.basic_property_pass].iterrows():
  mol=Chem.MolFromSmiles(r.standardized_smiles);matches=catalog.GetMatches(mol)
  alerts=[m.GetDescription() for m in matches]
  evaluated.append({'library_id':r.library_id,'structural_alert_evaluated':True,'structural_alerts':'|'.join(alerts),'structural_alert_pass':len(alerts)==0})
  if not alerts:
   entry=r.to_dict();entry['qed']=QED.qed(mol);entry['qualified_rank']=len(chosen)+1;chosen.append(entry)
  if len(chosen)==1000:break
 top=pd.DataFrame(chosen);write(top,'official_qualified_top1000.csv');write(pd.DataFrame(evaluated),'official_alert_evaluations.csv')
 admet=read(H/'02_admet_top300.csv');dock=read(H/'03_docking_top60.csv');online=read(H/'05_online_evidence.csv')

 for evidence_name, evidence in [('ADMET',admet),('docking',dock)]:
  structure_col='isomeric_smiles'
  identity=predict(evidence[structure_col].tolist())
  if not identity.status.eq('OK').all():raise RuntimeError(evidence_name+' has invalid structure')
  if not np.array_equal(identity.connectivity_key.to_numpy(),evidence.connectivity_key.to_numpy()):raise RuntimeError(evidence_name+' connectivity mismatch')
  if evidence_name=='ADMET':
   library_identity=scores.set_index('library_id').full_inchikey
   expected=evidence.library_id.map(library_identity)
   if not np.array_equal(expected.to_numpy(),identity.full_inchikey.to_numpy()):raise RuntimeError('ADMET full structure identity mismatch')
 # Keep all 60 evidence rows for an explicit coverage audit, not only favorable entries.
 diagnostic=predict(dock.isomeric_smiles.tolist());assert diagnostic.status.eq('OK').all()
 assert np.array_equal(diagnostic.connectivity_key.to_numpy(),dock.connectivity_key.to_numpy())
 assert np.array_equal(dock.library_id.map(scores.set_index('library_id').full_inchikey).to_numpy(),diagnostic.full_inchikey.to_numpy()), 'Docking full identity differs from library'
 merged=dock.copy();merged['historical_model_score']=merged.model_score
 merged=merged.drop(columns=['model_score'])
 for c in ['ranking_score','full_inchikey','max_training_similarity','descriptors_in_training_range','applicability_domain','fixed_threshold_class']:merged[c]=diagnostic[c].to_numpy()
 merged['model_version']=lock['version'];merged['in_current_qualified_top1000']=merged.library_id.isin(set(top.library_id))
 merged['current_qualified_rank']=merged.library_id.map(top.set_index('library_id').qualified_rank)
 # Complete InChIKey in public identity evidence must match the newly standardized structure.
 fields=['library_id','inchikey','product_identity_confirmed','availability_component','formal_product_id','supplier','product_url','supplier_status','known_GPR84_antagonist','known_GPR84_inactive','supplier_checked_utc','pubchem_cid','iupac_name','pubchem_title']
 oi=online[[c for c in fields if c in online]].copy()
 for c in oi.columns:
  if c!='library_id' and c in merged:merged=merged.drop(columns=c)
 merged=merged.merge(oi,on='library_id',how='left',validate='one_to_one')
 merged['exact_identity_match']=merged.full_inchikey.eq(merged.inchikey)
 train=read(ROOT/'data/training/training_structures.csv');g=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
 refs=[g.GetFingerprint(Chem.MolFromSmiles(s)) for s in train.loc[train.label.eq(1),'isomeric_smiles']]
 positive_scaffolds=set(train.loc[train.label.eq(1),'murcko_scaffold'])
 fps=[g.GetFingerprint(Chem.MolFromSmiles(s)) for s in diagnostic.standardized_smiles]
 merged['nearest_official_positive_similarity']=[max(DataStructs.BulkTanimotoSimilarity(fp,refs)) for fp in fps]
 merged['novel_to_official_positive_scaffold']=~merged.murcko_scaffold.isin(positive_scaffolds)
 # Same saved component formula, updated GPR84 AI scores and novelty reference only.
 merged['official_ai_component']=merged.ranking_score.rank(pct=True)
 merged['official_docking_component']=.7*merged.vina_best_kcal_mol.rank(ascending=False,pct=True)+.2*(1-merged.vina_seed_sd.clip(0,1.5)/1.5)+.1*merged.key_residue_interactions.fillna('').str.split(';').map(lambda v:min(len([x for x in v if x]),4)/4)
 merged['official_admet_component']=.3*merged.bbb_probability+.3*merged.adme_component+.3*merged.tox_safety_component+.1*merged.qed
 merged['official_novelty_component']=.6*(1-(merged.nearest_official_positive_similarity-.45).abs()/.55).clip(0,1)+.4*merged.novel_to_official_positive_scaffold.astype(float)
 merged['official_composite']=.30*merged.official_ai_component+.25*merged.official_docking_component+.25*merged.official_admet_component+.10*merged.official_novelty_component+.10*merged.availability_component
 merged['official_hard_gate']=merged.in_current_qualified_top1000 & merged.exact_identity_match & merged.product_identity_confirmed.fillna(False) & merged.mandatory_admet_complete & merged.pains_pass & merged.reactive_alert_pass & merged.tox_safety_component.ge(.35)&merged.adme_component.ge(.30)&merged.qed.ge(.30)&merged.geometry_pass & merged.selected_state_seed_count.eq(3)&~merged.known_GPR84_antagonist.fillna(False)&~merged.known_GPR84_inactive.fillna(False)&merged.fixed_threshold_class.eq(1)
 merged['toxicology_alert_DILI']=merged.DILI.ge(.5);merged['toxicology_alert_hERG']=merged.hERG.ge(.5)
 ranked=merged[merged.official_hard_gate].sort_values(['official_composite','library_id'],ascending=[False,True])
 selected=[];scaffolds=set();ood=0
 fp_map=dict(zip(merged.library_id,fps))
 for _,r in ranked.iterrows():
  outside=not r.applicability_domain
  if r.murcko_scaffold in scaffolds or (outside and ood>=1):continue
  if any(DataStructs.TanimotoSimilarity(fp_map[r.library_id],fp_map[v['library_id']])>=.55 for v in selected):continue
  entry=r.to_dict();entry['candidate_code']=f'GPR84-F{len(selected)+1:02}';entry['candidate_status']='COMPUTATIONAL_EVIDENCE_INTEGRATED_CANDIDATE';selected.append(entry);scaffolds.add(r.murcko_scaffold);ood+=int(outside)
  if len(selected)==5:break
 result=pd.DataFrame(selected);write(merged,'official_scored_existing_docking_pool60.csv');write(ranked,'official_evidence_covered_eligible_pool.csv');write(result,'official_provisional_candidates.csv')
 missing=top[~top.library_id.isin(set(admet.library_id))].copy();missing['pending']='ADMET predictions; then candidate-specific docking and exact product evidence';write(missing,'current_top1000_missing_ADMET.csv')
 coverage=top[['library_id','qualified_rank','ranking_score','applicability_domain']].copy();coverage['ADMET_available']=coverage.library_id.isin(set(admet.library_id));coverage['docking_available']=coverage.library_id.isin(set(dock.library_id));write(coverage,'current_top1000_evidence_coverage.csv')
 summary={'model_version':lock['version'],'qualified_top1000':len(top),'basic_property_pass':int(ordered.basic_property_pass.sum()),'structural_alert_evaluated':len(evaluated),'ADMET_available_in_current_top1000':int(coverage.ADMET_available.sum()),'docking_available_in_current_top1000':int(coverage.docking_available.sum()),'ADMET_missing_current_top1000':len(missing),'ADMET_available_in_current_top300':int(coverage.head(300).ADMET_available.sum()),'full_gate_evidence_covered':len(ranked),'provisional_candidate_count':len(result),'candidate_codes':result.candidate_code.tolist(),'candidate_library_ids':result.library_id.tolist(),'candidate_scope':'Complete-evidence subset intersecting the official qualified Top1000','weights':{'AI':.30,'docking':.25,'ADMET_BBB':.25,'novelty':.10,'availability':.10},'selection_constraints':{'unique_Murcko':True,'pairwise_Tanimoto_strictly_less_than':.55,'outside_official_AD_at_most':1},'old_scores_used_for_selection':False,'new_ADMET_or_docking_executed':False,'human_signoff':False,'quotation_or_purchase':False,'wet_experiment_scope':'computational screening, evidence integration and result verification','global_final_top5_confirmed':False}
 (S/'SCREENING_INTEGRATION_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False));verify_lock()
