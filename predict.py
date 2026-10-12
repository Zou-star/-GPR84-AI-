"""One command: regenerate official 100k scores and the five submission candidates."""
import argparse,sys
from pathlib import Path
sys.dont_write_bytecode=True
from src.util import ROOT
from src.workflow import run
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results')
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    if a.workers<1:p.error('--workers must be >= 1')
    run(a.output,a.workers)
if __name__=='__main__':main()
