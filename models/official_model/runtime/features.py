from .common import *
import rdkit
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator,MACCSkeys
def features(d):
    mols=[Chem.MolFromSmiles(s) for s in d.isomeric_smiles]
    def array(fps,size):
        a=np.zeros((len(fps),size),dtype=np.float32)
        for i,fp in enumerate(fps):DataStructs.ConvertToNumpyArray(fp,a[i])
        return a
    g4=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
    g6=rdFingerprintGenerator.GetMorganGenerator(radius=3,fpSize=2048,includeChirality=True)
    old=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
    ap=rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048,includeChirality=True)
    e4=array([g4.GetCountFingerprint(m) for m in mols],2048);e6=array([g6.GetCountFingerprint(m) for m in mols],2048)
    mac=array([MACCSkeys.GenMACCSKeys(m) for m in mols],167)
    desc=np.array([[standardize(s)[n] for n in DESC] for s in d.isomeric_smiles],dtype=np.float32)
    b4=array([g4.GetFingerprint(m) for m in mols],2048);b6=array([g6.GetFingerprint(m) for m in mols],2048)
    return {'F1':np.hstack([array([old.GetFingerprint(m) for m in mols],2048),desc]),'F2':np.hstack([e4,desc]),'F3':np.hstack([e6,desc]),'F4':np.hstack([mac,desc]),'F5':np.hstack([array([ap.GetCountFingerprint(m) for m in mols],2048),desc]),'F6':np.hstack([e4,e6,mac,desc]),'ECFP4':b4,'ECFP6':b6,'MACCS':mac}
def tanimoto(a,b):
    dot=a@b.T;den=a.sum(axis=1)[:,None]+b.sum(axis=1)[None,:]-dot
    return np.divide(dot,den,out=np.zeros_like(dot),where=den>0)
def main():
    files=[OUT/'02_splits/development.csv']+[p for p in (OUT/'02_splits').glob('window_*.csv')]
    for p in files:
        d=pd.read_csv(p)
        if not len(d):continue
        np.savez_compressed(OUT/'03_features'/f'{p.stem}.npz',**features(d),y=d.label.to_numpy(dtype=int))
        print(p.stem,len(d),flush=True)
    dump(OUT/'03_features/feature_schema.json',{'rdkit':rdkit.__version__,'descriptors':DESC,'F1':'nonchiral ECFP4 bits 2048 + 10 descriptors','F2':'chiral ECFP4 counts 2048 + 10 descriptors','F3':'chiral ECFP6 counts 2048 + 10 descriptors','F4':'MACCS 167 bits + 10 descriptors','F5':'chiral hashed AtomPair counts 2048 + 10 descriptors','F6':'chiral ECFP4 counts + ECFP6 counts + MACCS + 10 descriptors','Tanimoto':'chiral binary ECFP4/ECFP6 or MACCS without descriptors; pairwise intersection / union','standardization':'FragmentParent and Uncharger; stereochemistry preserved unless primary source explicitly reports racemic mixture'})
if __name__=='__main__':main()
