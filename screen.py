"""Arbitrary single or batch SMILES prediction; invalid rows are retained."""
import argparse,sys
from pathlib import Path
import pandas as pd
sys.dont_write_bytecode=True
from src.util import ROOT,output_file,write_csv,verify_inputs
from src.official import predict
def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--input',type=Path);g.add_argument('--smiles',nargs='+')
    p.add_argument('--output',type=Path,default=ROOT/'results/custom_scores.csv')
    a=p.parse_args();verify_inputs()
    if a.output.suffix.lower()!='.csv':p.error('--output must be a .csv file')
    try: destination=output_file(a.output)
    except ValueError as e:p.error(str(e))
    if a.input and a.input.resolve()==destination:p.error('Cannot overwrite the input CSV')
    if a.input:
        d=pd.read_csv(a.input,keep_default_na=False,skip_blank_lines=False,encoding='utf-8-sig')
        col=next((c for c in ['SMILES','smiles','isomeric_smiles'] if c in d),None)
        if col is None:p.error('CSV needs SMILES, smiles or isomeric_smiles column')
        values=d[col].tolist()
    else:values=a.smiles
    # Chunking bounds memory used by 400-reference applicability comparisons.
    parts=[predict(values[i:i+500]) for i in range(0,len(values),500)]
    result=pd.concat(parts,ignore_index=True) if parts else predict([])
    result['row_id']=range(len(result))
    write_csv(result,destination);print(destination.name)
if __name__=='__main__':main()
