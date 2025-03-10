import numpy as np
from rdkit import Chem
pt = Chem.GetPeriodicTable()

def is_molecule_confortable(atoms, coords, threshold_ratio=0.8):
    try:
        dist_mat = np.array([np.linalg.norm(coords - coord,axis=1) for coord in coords])
        for i in range(len(dist_mat)):
            dist_mat[i][i] = 100000

        atom_covR = np.array([pt.GetRcovalent(atom) for atom in atoms])
        aa_cov_dist_mat = np.array([[atomi+atomj for atomi in atom_covR] for atomj in atom_covR])

        return np.min(dist_mat - aa_cov_dist_mat * threshold_ratio) > 0
    except Exception as e:
        print(e)
        return False