import os, sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

BASE = Path("/project/home/p201228/gpcr_screening")

LIGANDS = {
    "common": [
        ("L01","O=c1[nH]c(OCCCOc2ccccc2Cl)nc2ncccc12"),
        ("L02","O=c1cc(CCCC2CC2)c2c(=O)[nH]c(C(F)F)nc2o1"),
        ("L03","CC1(CCCc2cc(=O)oc3nc(C(F)F)[nH]c(=O)c23)CC1"),
        ("L04","CCCCc1cc(=O)oc2[nH]c(=O)[nH]c(=O)c12"),
        ("L05","CC(C)CCc1cc(=O)oc2[nH]c(=O)[nH]c(=O)c12"),
        ("L06","O=c1[nH]c(=O)c2c(CCCC(F)F)cc(=O)oc2[nH]1"),
        ("L07","O=c1[nH]c(=O)c2c(CCC3CCC3)cc(=O)oc2[nH]1"),
        ("L08","O=c1[nH]c(=O)c2c(CCCC3CC3)cc(=O)oc2[nH]1"),
        ("L09","CC1(CCCc2cc(=O)oc3[nH]c(=O)[nH]c(=O)c23)CC1"),
        ("L10","CCc1cccc(C2(C)OC(C(=O)O)=CC2=O)c1"),
        ("L11","CC1(c2cccc(Br)c2)OC(C(=O)O)=CC1=O"),
        ("L12","CC1(c2cccc(C(F)(F)F)c2)OC(C(=O)O)=CC1=O"),
        ("L13","CC1(c2csc(Cl)c2)OC(C(=O)O)=CC1=O"),
        ("L14","CC1(c2cc(Br)cs2)OC(C(=O)O)=CC1=O"),
        ("L15","CC1(c2cccc(I)c2)OC(C(=O)O)=CC1=O"),
        ("L16","CC1(c2ccsc2)OC(C(=O)O)=CC1=O"),
        ("L17","Cc1ccc(C2(C)OC(C(=O)O)=CC2=O)s1"),
        ("L18","Cc1csc(C2(C)OC(C(=O)O)=CC2=O)c1"),
        ("L19","CC1(c2cc(F)cc(F)c2)OC(C(=O)O)=CC1=O"),
        ("L20","c1c(-c2nnn[nH]2)n[nH]c1C1CC1"),
        ("L21","CC(C)CNc1ccc(C(=O)O)cn1"),
        ("L22","CC(C)c1cc(C(=O)O)n[nH]1"),
        ("L23","CCCCCc1cc(C(=O)O)[nH]n1"),
        ("L24","CC1(c2ccccc2)OC(C(=O)O)=CC1=O"),
        ("L25","Cc1sc(C2(C)OC(C(=O)O)=CC2=O)cc1Br"),
        ("L26","CC1(c2cccc(Cl)c2)OC(C(=O)O)=CC1=O"),
        ("L27","CC1(c2csc(Br)c2)OC(C(=O)O)=CC1=O"),
        ("L28","CC1(c2c(F)cccc2F)OC(C(=O)O)=CC1=O"),
        ("L29","CC1(C2=CCCC2)OC(C(=O)O)=CC1=O"),
        ("L30","CC1(c2ccc(F)cc2)OC(C(=O)O)=CC1=O"),
        ("L31","CC1(c2ccc(F)c(F)c2)OC(C(=O)O)=CC1=O"),
        ("L32","CC(C)n1nnc2cc(C(=O)O)ccc21"),
        ("L33","CC1(c2ccc(F)cc2F)OC(C(=O)O)=CC1=O"),
        ("L34","CCCc1cc(-c2nnn[nH]2)n[nH]1"),
        ("L35","CCCCc1cc(-c2nnn[nH]2)n[nH]1"),
        ("L36","CCC(CC)Nc1ccc(C(=O)O)cn1"),
        ("L37","CCc1[nH]nc(-c2nnn[nH]2)c1F"),
        ("L38","CC1(c2ccc(Cl)s2)OC(C(=O)O)=CC1=O"),
        ("L39","Cc1cc(C2(C)OC(C(=O)O)=CC2=O)cs1"),
        ("L40","CCCCc1cc(C(=O)O)n[nH]1"),
        ("L41","CCC1(c2ccccc2)OC(C(=O)O)=CC1=O"),
        ("L42","O=C(O)c1cc(N(Cc2ccsc2)Cc2ccsc2)[nH]n1"),
        ("L43","O=C(O)c1cccnc1"),
    ]
}

def make_pdbqt(name, smiles, out_path):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol)

        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])

        if not ok:
            raise ValueError(err)

        out_path.write_text(pdbqt_str)
        print(f"[OK] {name}")
        return True

    except Exception as e:
        print(f"[ERROR] {name}: {e}")
        return False


out_dir = BASE / "ligands_common_43/pdbqt"
out_dir.mkdir(parents=True, exist_ok=True)

ok = fail = 0
for name, smi in LIGANDS["common_43"]:
    if make_pdbqt(name, smi, out_dir / f"{name}.pdbqt"):
        ok += 1
    else:
        fail += 1

print(f"\n{ok}/43 OK | {fail} fallos")
