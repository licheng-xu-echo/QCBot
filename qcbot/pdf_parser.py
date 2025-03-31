import fitz,re,json,os
from rdkit import Chem
import numpy as np
from collections import Counter
pt = Chem.GetPeriodicTable()

# API_KEY = "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" # from openrouter.ai
# API_KEY = "sk-U6CCn1011i1Fa2Ap979TLUq5qQdRWvEGK6DW25uppoBngsjU" # from tencent
def drop_page_id_line(page_lines,page_id,page_id_len=15,bias=3):
    while True:
        triggle = False
        for idx in [0,1]+[-1]:
            if idx >= len(page_lines) or -idx >= len(page_lines):
                continue
            match_flag = False
            for b in range(bias):
                if str(page_id-b) in page_lines[idx] and len(page_lines[idx]) <= page_id_len:
                    match_flag = True
                    break
                if str(page_id+b) in page_lines[idx] and len(page_lines[idx]) <= page_id_len:
                    match_flag = True
                    break
                if check_s_number_pattern(page_lines[idx]) and len(page_lines[idx]) <= page_id_len and not 'ts' in page_lines[idx].lower():
                    match_flag = True
                    break
            if match_flag:
                page_lines.pop(idx)
                triggle = True
                break
            if "supporting information" in page_lines[idx].lower():
                page_lines.pop(idx)
                triggle = True
                break
        if not triggle:
            break
    return page_lines

def drop_page_with_specific_elements(page_lines,elements):
    while True:
        triggle = False
        for idx in [0,1,2,3]+[-1,-2,-3,-4]:
            if idx >= len(page_lines) or -idx >= len(page_lines):
                continue
            for element in elements:
                if element in page_lines[idx]:
                    page_lines.pop(idx)
                    triggle = True
                    break
            if triggle:
                break
        if not triggle:
            break
    return page_lines

def clean_empty_lines(lines):
    new_lines = []
    for line in lines:
        if list(set(line)) == [""] or list(set(line)) == [" "] or line == "" or line == " ":
            continue
        new_lines.append(line)
    return new_lines

def drop_space(lines):
    new_lines = []
    for line in lines:
        new_lines.append(" ".join(line.split()))
    return new_lines
    
def check_s_number_pattern(s):
    # 正则表达式模式：匹配"S"后紧跟至少一个数字
    pattern = r'S\d+'
    return bool(re.search(pattern, s))

def extract_json_dict1(input_str):
    """
    从包含JSON的字符串中提取字典结构
    支持处理代码块包裹、转义字符和格式修正
    """
    try:
        # 清理代码块标记和转义字符（文献[2][6]）
        cleaned_str = re.sub(r'^```json\n|\n```$', '', input_str, flags=re.MULTILINE)
        cleaned_str = cleaned_str.replace('\\n', '\n').strip()
        
        # 处理可能存在的单引号问题（文献[5][7]）
        cleaned_str = cleaned_str.replace("'", '"')
        
        # 解析JSON（文献[1][4][8]）
        return json.loads(cleaned_str)
    
    except json.JSONDecodeError as e:
        # 自动修复常见格式错误（文献[7]）
        try:
            # 尝试提取有效JSON部分（文献[2]）
            matches = re.findall(r'\{.*?\}', cleaned_str, flags=re.DOTALL)
            if matches:
                return json.loads(max(matches, key=len))
            raise
        except:
            raise ValueError(f"JSON解析失败：{str(e)}") from None

def extract_json_dict2(input_str):
    return eval(input_str)

def extract_json_dict3(input_str):
    """
    增强版JSON提取函数，支持处理含特殊字符的键名
    """
    try:
        # 清理代码块标记和转义字符
        cleaned_str = re.sub(r'^```json\s*|\s*```$', '', input_str, flags=re.MULTILINE)
        cleaned_str = cleaned_str.replace('\\n', '\n').strip()
        
        # 处理转义单引号（保留键名中的合法单引号）
        cleaned_str = cleaned_str.replace('\\\'', "'")
        
        # 解析JSON
        return json.loads(cleaned_str)
    
    except json.JSONDecodeError as e:
        try:
            # 自动修复常见格式错误
            repaired_str = re.sub(r',(\s*[]}])', r'\1', cleaned_str)  # 移除尾部逗号
            repaired_str = re.sub(r'[\x00-\x1F]', '', repaired_str)   # 移除控制字符
            
            return json.loads(repaired_str)
        except:
            raise ValueError(f"JSON解析失败：{str(e)}") from None
        
def extract_json_dict4(input_str):
    """
    增强版JSON解析函数，修复单引号和特殊符号问题
    """
    try:
        # 1. 清理代码块标记
        cleaned_str = re.sub(r'^```json\s*|\s*```$', '', input_str, flags=re.MULTILINE)
        
        # 2. 处理单引号（分步替换）
        cleaned_str = (
            cleaned_str
            .replace("\\'", "'")     # 恢复转义单引号
            .replace("'", '"')       # 替换单引号为双引号
            .replace('[\"', '["')    # 修复双引号转义
            .replace('\"]', '"]')    # 修复数组闭合
        )
        
        # 3. 处理科学计数法数字（如1e-5）
        cleaned_str = re.sub(r'(\d)([-+])(\d)', r'\1e\2\3', cleaned_str)
        
        # 4. 自动修复常见格式错误
        repaired_str = re.sub(r',(\s*[\]}])', r'\1', cleaned_str)  # 移除尾部逗号
        repaired_str = re.sub(r'[\x00-\x1F]', '', repaired_str)   # 删除控制字符
        
        # 5. 解析前验证
        return json.loads(repaired_str)
    
    except json.JSONDecodeError as e:
        # 精确错误定位
        error_line = repaired_str.split('\n')[e.lineno - 1]
        error_ptr = ' ' * (e.colno - 1) + '^'
        print(f"解析失败于第 {e.lineno} 行:\n{error_line}\n{error_ptr}")
        raise

def extract_json_dict(input_str):
    try:
        return extract_json_dict1(input_str)
    except ValueError:
        try:
            return extract_json_dict2(input_str)
        except:
            try:
                return extract_json_dict3(input_str)
            except:
                try:
                    return extract_json_dict4(input_str)
                except:
                    raise ValueError("JSON解析失败，请检查输入格式。")
            
def merge_fake_molecules(table):
    """
    将所有标号为 "XXXX" 的分子数据合并到之前最近的非 "XXXX" 标号中
    参数:
        table (list): 包含字典的列表，每个字典表示分子数据
    返回:
        list: 合并后的列表，所有 "XXXX" 数据已迁移到最近的合法标号
    """
    last_valid_key = None
    last_valid_entry = None
    
    for entry in table:
        # 遍历字典键的副本，避免修改冲突
        for key in list(entry.keys()):
            if key == "XXXX":
                if last_valid_entry and last_valid_key:
                    # 合并到最近的合法标号
                    last_valid_entry[last_valid_key].extend(entry[key])
                del entry[key]  # 移除已处理的 "XXXX" 键
            else:
                # 更新最近的合法标号及对应字典
                last_valid_key = key
                last_valid_entry = entry
                
    return table

def symbol_pos_to_xyz_file(symbols,positions,xyz_file,title=""):
    with open(xyz_file, 'w') as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{title}\n")
        for symbol, pos in zip(symbols, positions):
            f.write(f"{symbol:4s} {pos[0]:15f} {pos[1]:15f} {pos[2]:15f}\n")
            
def get_charge_and_mult(atoms, charge=0, is_atom=False):
    """
    计算分子体系的电荷和自旋多重度
    参数：
        mol_input: 分子SMILES表达式 或 原子符号列表(如 ['C','H','H','H','H'])
        charge: 预设电荷（默认0）
        is_atom: 是否为单原子体系（需特殊处理）
    返回：
        (charge, multiplicity)
    """


    # 计算总电子数（原子序数之和 - 电荷）
    if isinstance(atoms[0],int):
        total_e = sum(atoms) - charge
    else:
        total_e = sum(pt.GetAtomicNumber(symbol) for symbol in atoms) - charge
    multiplicity = 1  # 默认闭壳层

    # 自旋多重度判断逻辑
    if total_e % 2 == 1:  # 奇数电子
        multiplicity = 2
    elif is_atom:  # 单原子特殊处理[1,4](@ref)
        atomic_num = pt.GetAtomicNumber(atoms[0])
        valence_e = get_valence_electrons(atomic_num)
        multiplicity = 2 * (valence_e % 2) + 1 if valence_e < 14 else 3
    elif total_e == 0:  # 特殊情况处理
        raise ValueError("电子总数不能为0")

    return (charge, multiplicity)

def generate_g16_input_file(gjf_file,atoms,coords,charge,multiplicity,method,jobtype,cpu,memory,title="Generated by QuantumCalcAgent",ckpt=False):
    information = [
    f'%nproc={cpu}',
    f'%mem={memory}GB',
    f'# {jobtype} {method}',
    "",
    title,
    "",
    f"{charge} {multiplicity}",
    ]
    if ckpt:
        information = [f"%chk={os.path.basename(gjf_file)[:-4]}.chk"] + information
    for atom, coord in zip(atoms, coords):
        information.append(f"{atom:5s} {coord[0]:15f} {coord[1]:15f} {coord[2]:15f}")
        
    information += ["",""]
    with open(gjf_file, 'w') as f:
        f.write("\n".join(information))

def get_valence_electrons(atomic_num):
    """获取原子的价电子数（简化版）"""
    period = pt.GetPeriod(atomic_num)
    group = pt.GetGroup(atomic_num)
    
    if period == 1: return 1
    if period == 2: return min(group, 8)
    if period >=3:  # 过渡金属简化处理[4](@ref)
        return pt.GetNOuterElecs(atomic_num)
    return 0

def split_into_continuous_groups(nums):
    if not nums:
        return []
    
    groups = []
    current_group = [nums[0]]
    
    for num in nums[1:]:
        if num[0] == current_group[-1][0] + 1:
            current_group.append(num)
        else:
            groups.append(current_group)
            current_group = [num]
    groups.append(current_group)  # 添加最后一个分组
    
    return groups

def find_header_and_footer_elements(nested_list, threshold=10):
    """
    在嵌套列表中找到出现次数超过指定阈值的元素
    
    参数:
    nested_list (list): 包含多个子列表的嵌套列表
    threshold (int): 频率阈值，默认10次
    
    返回:
    list: 超过阈值的元素列表
    """
    # 创建计数器对象
    element_counter = Counter()
    
    # 遍历每个子列表进行计数
    for sublist in nested_list:
        element_counter.update(sublist)
    
    # 筛选符合条件的元素
    return [element for element, count in element_counter.items() if count > threshold and len(element) > 5]

def parse_atom_coordinate(line):
    """
    解析原子坐标行，提取元素符号和三维坐标
    
    参数：
    line (str): 原子坐标字符串，格式示例：
                "Rh-0.84075200 -1.17783800 -0.44030400"
                "Rh -0.84075200 -1.17783800 -0.44030400"
                "Rh-0.84075200 -1.17783800 -0."
    
    返回：
    tuple: (元素符号, [x坐标, y坐标, z坐标])
           或 None（当格式不匹配时）
    """
    # 优化后的正则表达式模式
    pattern = re.compile(
        r'^\s*'                          # 行首空格
        r'([A-Za-z]{1,3})'               # 元素符号（组1）
        r'\d*'                           # 可选原子编号
        r'([-+]?\s*)?'                   # 符号与元素间的连接符（组2）
        r'([+-]?\d+\.?\d*)'              # X坐标（组3）
        r'\s+([+-]?\d+\.?\d*)'           # Y坐标（组4）
        r'\s+([+-]?\d+\.?\d*)'           # Z坐标（组5）
        r'\s*$'                          # 行尾空格
    )
    
    match = pattern.match(line)
    if not match:
        return None
    
    # 提取元素符号
    element = match.group(1)
    
    # 处理符号连接的特殊情况（如 Rh-0.5 → 需要合并符号和数值）
    x_coord = match.group(3)
    if match.group(2):  # 如果存在连接符号
        x_coord = match.group(2).strip() + x_coord
    
    return element, [x_coord, match.group(4), match.group(5)]
