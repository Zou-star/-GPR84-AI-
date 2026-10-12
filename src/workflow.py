"""Frozen official inference -> documented filters -> evidence-based ranking."""
import time, platform
import numpy as np,pandas as pd
from rdkit import Chem
from src.util import ROOT,OFF,sha,read_json,dump,write_csv,output_dir,verify_inputs,scoring_signature
from src.official import predict,verify_lock
def export_results(dest, demo=False):
    cfg=read_json(ROOT/'config.json')
    d=pd.read_csv(dest/'official_provisional_candidates.csv',float_precision='round_trip')
    result=[]
    for rank,(_,r) in enumerate(d.iterrows(),1):
        pose='data/evidence/structures/docking/candidates/'+r.library_id+'/best_pose.sdf'
        # Open through Python so RDKit also works with Windows Chinese paths.
        with (ROOT/pose).open('rb') as stream:
            mol=next(iter(Chem.ForwardSDMolSupplier(stream,removeHs=True)))
        if mol is None or Chem.MolToInchiKey(mol).split('-')[0]!=r.connectivity_key:
            raise RuntimeError('Pose/parent connectivity mismatch: '+r.library_id)
        result.append(dict(candidate_id=r.candidate_code,rank=rank,track=cfg['track'],target='human GPR84',
            library_id=r.library_id,SMILES=r.isomeric_smiles,full_inchikey=r.full_inchikey,
            ranking_score=float(r.ranking_score),fixed_threshold=float(verify_lock()['threshold']),
            fixed_threshold_class=int(r.fixed_threshold_class),applicability_domain=bool(r.applicability_domain),
            max_training_similarity=float(r.max_training_similarity),composite_score=float(r.official_composite),
            BBB_probability=float(r.bbb_probability),ADMET_component=float(r.official_admet_component),
            DILI=float(r.DILI),hERG=float(r.hERG),AMES=float(r.AMES),
            Vina_best_kcal_mol=float(r.vina_best_kcal_mol),Vina_seed_sd=float(r.vina_seed_sd),
            structure_file=pose,receptor_file='data/evidence/structures/docking/9WIJ.cif',
            supplier=r.supplier,product_id=r.formal_product_id,product_url=r.product_url,
            model_version=cfg['model_version'],run_version=cfg['run_version'],
            candidate_status=r.candidate_status,independent_acceptance=False,
            notes='基于正式模型、ADMET/BBB预测、三种子对接和结构身份核对的计算筛选候选；证据范围见Model Card'))
    write_csv(pd.DataFrame(result),dest/'results.csv')
    return [r['candidate_id'] for r in result]
def run(dest, workers=4, demo=False):
    start=time.perf_counter();dest=output_dir(dest);count=verify_inputs();lock=verify_lock()
    if demo:
        if dest==ROOT/'results':raise ValueError('Demo must use a separate output directory')
        # The 60-row demo retains a frozen Top1000 membership table so that it
        # demonstrates the same decision rules. It is not a fresh 100k screen.
        input_df=pd.read_csv(ROOT/'data/evidence/03_docking_top60.csv',float_precision='round_trip')
        pr=predict(input_df.isomeric_smiles.tolist())
        if not pr.status.eq('OK').all():raise RuntimeError('Invalid demonstration structure')
        for c in ['ranking_score','full_inchikey','connectivity_key','applicability_domain']:
            input_df[c]=pr[c].to_numpy()
        top=pd.read_csv(ROOT/'results/official_qualified_top1000.csv',float_precision='round_trip')
        # Demo uses source-frozen qualifying membership rather than pretending
        # that 60 rows reproduce the global Top1000 search.
        from src import evidence_ranking as ranking
        ranking.S=dest
        admet=pd.read_csv(ROOT/'data/evidence/02_admet_top300.csv',float_precision='round_trip')
        dock=pd.read_csv(ROOT/'data/evidence/03_docking_top60.csv',float_precision='round_trip')
        if not np.array_equal(pr.full_inchikey.to_numpy(),input_df.full_inchikey.to_numpy()):raise RuntimeError('Demo identity mismatch')
        # A common integration routine is used below with a reconstructed
        # universe of source rows and frozen Top1000 membership flags.
        saved=pd.read_parquet(ROOT/'results/official_100k_scores.parquet')
        baseline=read_json(ROOT/'data/reference/scoring_baseline.json')
        signature=scoring_signature(saved)
        for key in signature:
            if signature[key]!=baseline[key]:raise RuntimeError('Saved demo global membership baseline differs: '+key)
        needed=set(top.library_id)|set(admet.library_id)|set(dock.library_id)
        subset=saved[saved.library_id.isin(needed)].copy()
        overlap=subset.library_id.isin(set(input_df.library_id))
        subset.loc[overlap,'ranking_score']=subset.loc[overlap,'library_id'].map(input_df.set_index('library_id').ranking_score)
        subset.to_parquet(dest/'official_100k_scores.parquet',index=False)
        # Evidence ranking rebuilds qualifying membership from this frozen
        # eligible Top1000 subset; excluded PAINS rows are unnecessary here.
        ranking.main(dest)
    else:
        from src.library_scoring import main as score
        from src.evidence_ranking import main as integrate
        score(dest,workers);integrate(dest)
    ids=export_results(dest,demo)
    summary=read_json(dest/'SCREENING_INTEGRATION_SUMMARY.json')
    dump(dest/'run_receipt.json',{'status':'COMPLETE','mode':'small_demo_with_frozen_global_membership' if demo else 'fresh_official_100k_scoring_with_frozen_downstream',
        'model_version':lock['version'],'model_sha256':sha(OFF/'model.joblib'),
        'library_reference_sha256':read_json(ROOT/'data/library/library_manifest.json')['original_snapshot_sha256'],
        'immutable_files_verified':count,'candidate_ids':ids,'elapsed_seconds':time.perf_counter()-start,
        'environment':{'python':platform.python_version(),'platform':platform.platform()},
        'model_refitted':False,'new_ADMET_or_docking_executed':False,'global_final_top5_confirmed':False,
        'independent_acceptance':False,'evidence_summary':summary,
        'outputs_sha256':{p.name:sha(p) for p in dest.iterdir() if p.is_file() and p.name!='run_receipt.json'}})
    print('RESULTS',(dest/'results.csv').relative_to(ROOT) if dest.is_relative_to(ROOT) else dest.name+'/results.csv',flush=True);return dest
