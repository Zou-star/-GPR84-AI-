"""Small, offline demonstration using frozen full-library qualification membership."""
import argparse,sys
from pathlib import Path
sys.dont_write_bytecode=True
from src.util import ROOT
from src.workflow import run
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'results/demo')
    a=p.parse_args();run(a.output,demo=True)
if __name__=='__main__':main()
