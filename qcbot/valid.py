import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

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

def check_lonely_atoms(atom_types, coords, atom_type="H", threshold_dist=2.0):
    """
    检查分子中所有氢原子周围2埃内是否存在非氢原子
    
    参数:
    atom_types (list): 原子类型列表，如 ['H', 'C', 'O']
    coords (np.ndarray): 原子坐标数组，形状为(N,3)
    
    返回:
    bool: 如果所有氢原子周围都存在非氢原子返回True，否则返回False
    """
    coords = np.asarray(coords)
    h_indices = [i for i, atom in enumerate(atom_types) if atom == atom_type]
    
    # 如果没有氢原子直接返回False
    if not h_indices:
        return False
    
    atom_types_array = np.array(atom_types)
    
    for h_idx in h_indices:
        # 计算所有原子到当前氢原子的距离
        distances = np.linalg.norm(coords - coords[h_idx], axis=1)
        
        # 创建条件掩码：非氢原子、距离<=2埃、排除自身
        condition = (
            #(atom_types_array != 'H') & 
            (distances <= threshold_dist) & 
            (np.arange(len(atom_types)) != h_idx)
        )
        
        # 如果没有满足条件的原子，返回False
        if not np.any(condition):
            return False
    
    return True
def find_isolated_atoms(atoms,coords, tolerance=1.2):
    """
    检测XYZ文件中是否存在孤立原子
    返回：(是否存在问题, 孤立原子列表)
    """
    n_atoms = len(atoms)
    
    # 构建半径矩阵
    radius_matrix = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(n_atoms):
            r1 = pt.GetRcovalent(atoms[i])
            r2 = pt.GetRcovalent(atoms[j])
            radius_matrix[i,j] = (r1 + r2) * tolerance
    
    # 计算距离矩阵
    dist_matrix = cdist(coords, coords, 'euclidean')
    np.fill_diagonal(dist_matrix, np.inf)  # 忽略自身
    
    # 查找孤立原子
    isolated = []
    for i in range(n_atoms):
        if not np.any(dist_matrix[i] < radius_matrix[i]):
            isolated.append((i+1, atoms[i]))  # 原子编号从1开始
            
    return len(isolated) > 0, isolated