"""Fixed-parameter reproduction; no search and no replacement of official weights.

The default is a read-only training-input check. An explicit --reproduce trains
a separate comparison model. Training-set scores are diagnostics, not validation.
"""
import argparse,sys,time
from pathlib import Path
import joblib,numpy as np,pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from xgboost import XGBClassifier
sys.dont_write_bytecode=True
from src.util import ROOT,OFF,read_json,dump,sha,output_dir,verify_inputs,write_csv
from src.official import load_model
def check():
    verify_inputs();lock,bundle=load_model()
    d=pd.read_csv(ROOT/'data/training/training_structures.csv',float_precision='round_trip')
    m=pd.read_csv(ROOT/'data/training/training_members.csv');z=np.load(ROOT/'data/training/development.npz')
    if len(d)!=400 or int(d.label.sum())!=326 or int((d.label==0).sum())!=74:raise RuntimeError('Member counts changed')
    if not np.array_equal(d.compound_id,m.compound_id) or not np.array_equal(d.label,z['y']):raise RuntimeError('Feature/member order mismatch')
    if d.compound_id.duplicated().any() or d.full_inchikey.duplicated().any():raise RuntimeError('Training identity duplicate')
    # Distinct stereoisomers intentionally share connectivity keys. Keep all
    # frozen members and check grouping rather than deleting those structures.
    groups=pd.read_csv(ROOT/'data/training/publication_groups.csv')
    assert np.array_equal(groups.compound_id,d.compound_id)
    for _,indices in d.groupby('connectivity_key').groups.items():
        if groups.iloc[list(indices)].group.nunique()!=1:raise RuntimeError('Stereoisomers span publication groups')
    # Compare cached features against the unchanged frozen chemistry routines.
    sys.path.insert(0,str(OFF))
    from runtime.features import features
    fresh=features(d)['F4']
    if not np.array_equal(fresh,z['F4'],equal_nan=True):raise RuntimeError('Frozen training feature mismatch')
    if fresh.shape!=(400,177):raise RuntimeError('Unexpected feature dimensions')
    return lock,bundle,d,z
def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group()
    g.add_argument('--check-only',action='store_true');g.add_argument('--reproduce',action='store_true')
    p.add_argument('--output',type=Path,default=ROOT/'results/training_reproduction');a=p.parse_args()
    start=time.perf_counter();lock,bundle,d,z=check()
    if not a.reproduce:
        print('VERIFIED: 400 members; 326 positive; 74 negative; MACCS+10 = 177; no fitting');return
    dest=output_dir(a.output)
    if (dest/'reproduced_model.joblib').exists():raise FileExistsError('No overwrite of prior reproduction')
    params=dict(lock['parameters']);params.update(scale_pos_weight=74/326,n_jobs=2,random_state=42,tree_method='hist',eval_metric='logloss')
    model=make_pipeline(SimpleImputer(strategy='median'),XGBClassifier(**params))
    # Identical median imputation and fixed 500-round estimator. No early stopping
    # or validation labels are used for this final 400-member reproduction.
    model.fit(z['F4'],z['y'])
    reproduced=dict(bundle);reproduced['model']=model;reproduced['version']='NONOFFICIAL-FIXED-PARAMETER-REPRODUCTION'
    joblib.dump(reproduced,dest/'reproduced_model.joblib')
    old=bundle['model'].predict_proba(z['F4'])[:,1];new=model.predict_proba(z['F4'])[:,1]
    out=d[['compound_id','label']].copy();out['official_training_diagnostic_score']=old;out['reproduced_training_diagnostic_score']=new
    write_csv(out,dest/'training_score_comparison.csv')
    dump(dest/'reproduction_receipt.json',{'seed':42,'parameters':params,'N':400,'positive':326,'negative':74,
        'model_sha256':sha(dest/'reproduced_model.joblib'),'max_training_score_difference':float(np.max(np.abs(old-new))),
        'training_scores_bitwise_equal':bool(np.array_equal(old,new)),'elapsed_seconds':time.perf_counter()-start,
        'official_model_replaced':False,'validation_rerun':False,'scope':'fixed-parameter reproduction only; not new model selection or validation'})
    print(dest)
if __name__=='__main__':main()
