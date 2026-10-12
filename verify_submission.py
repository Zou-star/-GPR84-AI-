"""Verify immutable inputs, all-library numerical baseline and candidate output.

Verification is separate from candidate generation; --results-dir selects the
actual generated result directory, while --output selects the verification log.
"""
import argparse,sys,json,tempfile,shutil,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,matthews_corrcoef,balanced_accuracy_score
sys.dont_write_bytecode=True
from src.util import ROOT,OFF,verify_inputs,read_json,dump,sha,output_dir,output_file,scoring_signature
from src.official import predict,verify_lock,classify

def run(dest,results_dir,recompute=False):
    start=time.perf_counter();count=verify_inputs();lock=verify_lock()
    results_dir=Path(results_dir).resolve();dest=output_dir(dest)
    members=pd.read_csv(ROOT/'data/training/training_members.csv')
    assert len(members)==400 and int(members.label.sum())==326 and int((members.label==0).sum())==74
    assert sha(ROOT/'data/training/training_members.csv')==lock['training_member_hash']
    assert not members.full_inchikey.duplicated().any() and not members.compound_id.duplicated().any()
    groups=pd.read_csv(ROOT/'data/training/publication_groups.csv')
    assert np.array_equal(members.compound_id,groups.compound_id)
    for _,ix in members.groupby('connectivity_key').groups.items():
        assert groups.iloc[list(ix)].group.nunique()==1,'Connected stereoisomers must stay grouped'
    for fold in read_json(ROOT/'data/training/splits/folds.json')['publication']:
        tr,va=fold['train'],fold['validation']
        assert not(set(groups.iloc[tr].group)&set(groups.iloc[va].group))
        assert not(set(members.iloc[tr].connectivity_key)&set(members.iloc[va].connectivity_key))
    values=['CCO','invalid_SMILES','',None,'CCO.Cl','C[C@H](O)C(=O)O','C[C@@H](O)C(=O)O','CCO']
    out=predict(values)
    assert out.row_id.tolist()==list(range(8))
    assert out.status.tolist()==['OK','invalid_SMILES','invalid_SMILES','invalid_SMILES','OK','OK','OK','OK']
    assert out.iloc[0].ranking_score==out.iloc[7].ranking_score==out.iloc[4].ranking_score
    assert out.iloc[5].standardized_smiles!=out.iloc[6].standardized_smiles
    assert out.iloc[5].connectivity_key==out.iloc[6].connectivity_key
    assert classify(lock['threshold'],lock['threshold'])==1
    assert classify(np.nextafter(lock['threshold'],-np.inf),lock['threshold'])==0
    protected=['config.json','README.md','IMMUTABLE_SHA256.json','checksums.sha256',
               'models/official_model/model.joblib','data/library/public_purchasable_100k.csv',
               'src/util.py','logs/history/XGBoost_F4_trials.csv']
    for name in protected:
        try:output_file(ROOT/name)
        except ValueError:pass
        else:raise RuntimeError('Protected output accepted: '+name)
    with tempfile.TemporaryDirectory(prefix='tamper_',dir=str(dest)) as temporary:
        b=Path(temporary)/'official_model';shutil.copytree(OFF,b)
        with (b/'model.joblib').open('ab') as f:f.write(b'INTENTIONAL_TEST_CORRUPTION')
        try:verify_lock(b)
        except RuntimeError as e:assert 'hash mismatch' in str(e).lower()
        else:raise RuntimeError('Tampered model accepted')
    oof=pd.read_csv(ROOT/'data/training/development_publication_OOF.csv',float_precision='round_trip')
    y=oof.label.to_numpy();s=oof.score.to_numpy();p=s>=lock['threshold']
    estimates={'ROC_AUC':float(roc_auc_score(y,s)),'MCC':float(matthews_corrcoef(y,p)),
               'Balanced_Accuracy':float(balanced_accuracy_score(y,p))}
    for key,v in estimates.items():assert abs(v-lock['development_publication_OOF'][key])<1e-12
    r=pd.read_csv(results_dir/'results.csv',float_precision='round_trip')
    ref=pd.read_csv(ROOT/'data/reference/official_provisional_candidates.csv',float_precision='round_trip')
    cfg=read_json(ROOT/'config.json')
    assert len(r)==5 and r.library_id.tolist()==ref.library_id.tolist()
    assert r.candidate_id.tolist()==ref.candidate_code.tolist()
    assert r.track.eq(cfg['track']).all() and r.model_version.eq(lock['version']).all()
    assert r.run_version.eq(cfg['run_version']).all() and r.target.eq('human GPR84').all()
    assert np.array_equal(r.ranking_score.to_numpy(),ref.ranking_score.to_numpy())
    assert np.allclose(r.composite_score.to_numpy(),ref.official_composite.to_numpy(),rtol=0,atol=1e-12)
    assert not r.independent_acceptance.any()
    assert not any(c.startswith('wet_experiment') for c in r.columns)
    for relative in r.structure_file:
        target=(ROOT/relative).resolve()
        assert target.is_relative_to(ROOT) and target.is_file()
    scores=pd.read_parquet(results_dir/'official_100k_scores.parquet')
    signature=scoring_signature(scores);baseline=read_json(ROOT/'data/reference/scoring_baseline.json')
    for key in signature:assert signature[key]==baseline[key],'100k numerical baseline differs: '+key
    # Real CLI checks preserve invalid, blank and duplicate CSV records, including
    # the boundary between the two 500-row processing chunks.
    with tempfile.TemporaryDirectory(prefix='cli_',dir=str(dest)) as temporary:
        temp=Path(temporary)
        sample=pd.DataFrame({'SMILES':['CCO']*499+['invalid_SMILES','','CCO.Cl','CCO']})
        sample.to_csv(temp/'input.csv',index=False,encoding='utf-8-sig')
        cli=subprocess.run([sys.executable,str(ROOT/'screen.py'),'--input',str(temp/'input.csv'),
                            '--output',str(temp/'output.csv')],cwd=temp,text=True,encoding='utf-8',capture_output=True)
        if cli.returncode:raise RuntimeError('CLI inference failed: '+cli.stderr)
        actual=pd.read_csv(temp/'output.csv',float_precision='round_trip',keep_default_na=False)
        assert actual.row_id.tolist()==list(range(len(sample)))
        assert actual.status.tolist()==['OK']*499+['invalid_SMILES','invalid_SMILES','OK','OK']
        assert np.array_equal(actual.loc[actual.status.eq('OK'),'ranking_score'].astype(float),
                              np.full(501,out.iloc[0].ranking_score))
        before=sha(ROOT/'config.json')
        refusal=subprocess.run([sys.executable,str(ROOT/'screen.py'),'--smiles','CCO',
                                '--output',str(ROOT/'config.json')],cwd=temp,capture_output=True)
        assert refusal.returncode!=0 and sha(ROOT/'config.json')==before
    if recompute:
        from src.development import compute
        with tempfile.TemporaryDirectory(prefix='bootstrap_',dir=str(dest)) as temporary:
            compute(temporary)
            assert read_json(Path(temporary)/'development_metrics.json')==read_json(ROOT/'data/training/development_metrics.json')
            assert sha(Path(temporary)/'bootstrap_replicates.csv')==sha(ROOT/'data/training/bootstrap_replicates.csv')
    result={'status':'PASSED','model_version':lock['version'],'run_version':cfg['run_version'],
        'model_sha256':sha(OFF/'model.joblib'),'immutable_files_verified':count,
        'training_counts':{'N':400,'positive':326,'negative':74},
        'fresh_100k_scores_exactly_match_reference':True,'Top5_identity_order_and_values_match':True,
        'input_edge_cases_passed':True,'CLI_blank_invalid_duplicate_chunk_order_passed':True,
        'protected_output_paths_rejected':True,'threshold_equality_positive':True,
        'tampered_model_rejected':True,'development_recomputed_metrics':estimates,
        'full_2000x2_bootstrap_recomputed':recompute,'official_weights_changed':False,
        'results_location':'project/'+results_dir.relative_to(ROOT).as_posix() if results_dir.is_relative_to(ROOT) else 'external_results_directory',
        'elapsed_seconds':time.perf_counter()-start,'independent_acceptance':False}
    dump(dest/'verification.json',result);print(json.dumps(result,ensure_ascii=False));return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir',type=Path,default=ROOT/'results')
    p.add_argument('--output',type=Path,default=ROOT/'results/verification')
    p.add_argument('--recompute-development',action='store_true')
    a=p.parse_args();run(a.output,a.results_dir,a.recompute_development)
