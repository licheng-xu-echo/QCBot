import numpy as np
from collections import defaultdict
periodic_table = [
    # 周期 0（特殊处理）
    #[None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None],
    
    # 第1周期（2元素）
    ['H', None, None, None, None, None, None, None, None, None, None, None, None,None, None, None, None, 'He'],
    
    # 第2周期（8元素）
    ['Li', 'Be', None, None, None, None, None, None, None, None, None, None, 'B', 'C', 'N', 'O', 'F', 'Ne'],
    
    # 第3周期（8元素）
    ['Na', 'Mg', None, None, None, None, None, None, None, None, None, None, 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
    
    # 第4周期（含过渡金属）
    ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
    
    # 第5周期
    ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe'],
    
    # 第6周期（含镧系占位符）
    ['Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'],
    
    # 第7周期（含锕系占位符）
    ['Fr', 'Ra', 'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
]

# 镧系和锕系展开表（单独处理）
lanthanides = ['La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']
actinides = ['Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

def generate_metal_count_array(metal_lst):
    # 创建计数字典
    counter = defaultdict(int)
    for elem in metal_lst:
        counter[elem] += 1

    # 初始化二维数组
    arr = np.zeros((7, 18), dtype=int)
    
    # 主表填充
    for i, row in enumerate(periodic_table):
        for j, elem in enumerate(row):
            if elem in counter:
                arr[i][j] = counter.pop(elem)
    
    # 处理镧系/锕系
    lanthanide_count = sum(counter[e] for e in lanthanides if e in counter)
    actinide_count = sum(counter[e] for e in actinides if e in counter)
    
    # 镧系位置（第6周期第2列）
    arr[5][1] = lanthanide_count
    # 锕系位置（第7周期第2列）
    arr[6][1] = actinide_count

    return arr