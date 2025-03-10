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
    
def check_hydrogens(atom_types, coords, threshold_dist=2.0):
    """
    检查分子中所有氢原子周围2埃内是否存在非氢原子
    
    参数:
    atom_types (list): 原子类型列表，如 ['H', 'C', 'O']
    coords (np.ndarray): 原子坐标数组，形状为(N,3)
    
    返回:
    bool: 如果所有氢原子周围都存在非氢原子返回True，否则返回False
    """
    coords = np.asarray(coords)
    h_indices = [i for i, atom in enumerate(atom_types) if atom == 'H']
    
    # 如果没有氢原子直接返回False
    if not h_indices:
        return False
    
    atom_types_array = np.array(atom_types)
    
    for h_idx in h_indices:
        # 计算所有原子到当前氢原子的距离
        distances = np.linalg.norm(coords - coords[h_idx], axis=1)
        
        # 创建条件掩码：非氢原子、距离<=2埃、排除自身
        condition = (
            (atom_types_array != 'H') & 
            (distances <= threshold_dist) & 
            (np.arange(len(atom_types)) != h_idx)
        )
        
        # 如果没有满足条件的原子，返回False
        if not np.any(condition):
            return False
    
    return True