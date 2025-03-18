import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation
pt = Chem.GetPeriodicTable()

class MoleculePose:
    def __init__(self, atoms, coordinates, use_mass_weighted=False, reference_atom_idx=0):
        """
        输入参数：
        - atoms: list[str]，原子类型列表（如 ["C", "O", "H", ...]）
        - coordinates: np.array，形状为 (n_atoms, 3) 的原子坐标数组
        - use_mass_weighted: 是否使用原子质量加权计算质心和主方向（默认用几何中心）
        """
        self.atoms = atoms
        self.coordinates = coordinates
        self.use_mass_weighted = use_mass_weighted
        
        # 计算质心和主方向
        self.reference_atom_idx = reference_atom_idx  # 新增
        if len(self.atoms) > 1:
            self.center = self._compute_center()
            self.principal_axes = self._compute_principal_axes()
        else:
            self.center = self.coordinates[0]
            self.principal_axes = np.eye(3)  # 单原子分子，主方向为单位矩阵
    
    def _compute_center(self):
        """计算质心（几何中心或质量加权中心）"""
        if self.use_mass_weighted:
            masses = np.array([pt.GetAtomicWeight(atom) for atom in self.atoms])
            weighted_coords = self.coordinates * masses[:, np.newaxis]
            return np.sum(weighted_coords, axis=0) / np.sum(masses)
        else:
            return np.mean(self.coordinates, axis=0)
    
    def _compute_principal_axes(self):
        """通过 PCA 计算主方向，并修正符号"""
        centered = self.coordinates - self.center
        cov = np.cov(centered.T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        order = np.argsort(eigenvalues)[::-1]
        axes = eigenvectors[:, order]
        
        # ------------------------- 新增：统一主方向符号 -------------------------
        # 选择参考原子相对于质心的方向
        reference_vector = centered[self.reference_atom_idx]
        if np.linalg.norm(reference_vector) < 1e-10:
            reference_vector = np.array([1.0, 0.0, 0.0])  # 若参考原子在质心，使用默认方向
        
        # 对每个主方向，确保其与参考向量夹角小于 90 度
        for i in range(axes.shape[1]):
            axis = axes[:, i]
            # 计算当前主方向与参考向量的点积
            dot = np.dot(axis, reference_vector)
            # 如果点积为负，反转方向
            if dot < 0:
                axes[:, i] = -axis
        # ---------------------------------------------------------------------
        
        return axes
    
    #-----------------------------------------
    # 朝向表示方法
    #-----------------------------------------
    def get_center(self):
        """返回质心位置向量 [x, y, z]"""
        return self.center
    
    def get_principal_direction(self):
        """返回第一个主方向（单位向量）"""
        return self.principal_axes[:, 0]
    
    def get_rotation_matrix(self):
        """返回旋转矩阵（3x3），将分子对齐到主方向"""
        return self.principal_axes.T  # 列向量为基向量
    
    def get_quaternion(self):
        """返回四元数 [w, x, y, z]（从原始坐标系到主方向）"""
        rotation = Rotation.from_matrix(self.get_rotation_matrix())
        return rotation.as_quat()  # 格式 [x, y, z, w]
    
    def get_axis_angle(self):
        """返回轴角表示 [ax, ay, az, angle_radians]"""
        rotation = Rotation.from_matrix(self.get_rotation_matrix())
        return rotation.as_rotvec()  # 格式 [ax, ay, az]，角度为向量模长
    
    def get_euler_angles(self, sequence='zyx'):
        """返回欧拉角（弧度），默认顺序 Z-Y-X"""
        rotation = Rotation.from_matrix(self.get_rotation_matrix())
        return rotation.as_euler(sequence)
    
    def get_two_ortho_vectors(self):
        """返回前向和向上两个正交向量"""
        forward = self.principal_axes[:, 0]
        up = self.principal_axes[:, 1]
        return forward, up

    #-----------------------------------------
    # 可视化与验证
    #-----------------------------------------
    def align_molecule(self):
        """将分子坐标对齐到主方向，返回新坐标"""
        rotation_matrix = self.get_rotation_matrix()
        centered = self.coordinates - self.center
        return np.dot(centered, rotation_matrix.T)
    
    def restore_molecule(self, tgt_center, tgt_rot_mat):
        #!!!!!
        return np.dot(np.dot(self.coordinates-self.center, self.get_rotation_matrix().T), tgt_rot_mat) + tgt_center
    
    def random_transform(self, max_translation=5.0, max_rotation_deg=360.0, 
                        seed=None, rotate_around_center=True):
        """
        对分子进行随机旋转和平移
        
        参数：
            max_translation (float): 最大平移距离（单位：埃，默认5.0）
            max_rotation_deg (float): 最大旋转角度（单位：度，默认360.0）
            seed (int): 随机数种子（默认None）
            rotate_around_center (bool): 是否绕质心旋转（默认True）
        
        返回：
            np.ndarray: 变换后的坐标数组（形状与输入相同）
        """
        # 参数检查
        if max_translation < 0 or max_rotation_deg < 0:
            raise ValueError("平移和旋转参数必须为非负数")

        # 创建独立随机生成器（避免影响全局状态）
        rng = np.random.default_rng(seed)

        # 生成随机平移向量
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)  # 归一化
        translation = max_translation * direction

        # 生成随机旋转
        rotation = Rotation.from_rotvec(
            rng.uniform(-1, 1, 3) * np.radians(max_rotation_deg)
        )

        # 获取坐标副本（避免修改原始数据）
        new_coords = np.copy(self.coordinates)

        # 应用旋转
        if rotate_around_center:
            # 绕质心旋转
            centered = new_coords - self.center
            rotated = rotation.apply(centered)
            new_coords = rotated + self.center
        else:
            # 绕原点旋转
            new_coords = rotation.apply(new_coords)

        # 应用平移
        new_coords += translation

        return new_coords