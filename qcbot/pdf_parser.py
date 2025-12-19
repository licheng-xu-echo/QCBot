import fitz,re,json,os
from rdkit import Chem
import numpy as np
from collections import Counter
pt = Chem.GetPeriodicTable()

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
    # Regular expression pattern: Matches "S" followed by at least one digit
    pattern = r'S\d+'
    return bool(re.search(pattern, s))

def extract_json_dict1(input_str):
    """
    Extract the dictionary structure from a string containing JSON
    Support handling code block wrapping, escape characters, and formatting correction
    """
    try:

        cleaned_str = re.sub(r'^```json\n|\n```$', '', input_str, flags=re.MULTILINE)
        cleaned_str = cleaned_str.replace('\\n', '\n').strip()
        cleaned_str = cleaned_str.replace("'", '"')
        return json.loads(cleaned_str)
    
    except json.JSONDecodeError as e:
        
        try:
            matches = re.findall(r'\{.*?\}', cleaned_str, flags=re.DOTALL)
            if matches:
                return json.loads(max(matches, key=len))
            raise
        except:
            raise ValueError(f"JSON parse failed: {str(e)}") from None

def extract_json_dict2(input_str):
    return eval(input_str)

def extract_json_dict3(input_str):

    try:
        cleaned_str = re.sub(r'^```json\s*|\s*```$', '', input_str, flags=re.MULTILINE)
        cleaned_str = cleaned_str.replace('\\n', '\n').strip()
        
        cleaned_str = cleaned_str.replace('\\\'', "'")
        
        return json.loads(cleaned_str)
    
    except json.JSONDecodeError as e:
        try:
            
            repaired_str = re.sub(r',(\s*[]}])', r'\1', cleaned_str) 
            repaired_str = re.sub(r'[\x00-\x1F]', '', repaired_str)   
            
            return json.loads(repaired_str)
        except:
            raise ValueError(f"JSON parse failed: {str(e)}") from None
        
def extract_json_dict4(input_str):

    try:

        cleaned_str = re.sub(r'^```json\s*|\s*```$', '', input_str, flags=re.MULTILINE)
        
        cleaned_str = (
            cleaned_str
            .replace("\\'", "'")     
            .replace("'", '"')       
            .replace('[\"', '["')    
            .replace('\"]', '"]')    
        )
        

        cleaned_str = re.sub(r'(\d)([-+])(\d)', r'\1e\2\3', cleaned_str)
        

        repaired_str = re.sub(r',(\s*[\]}])', r'\1', cleaned_str)  
        repaired_str = re.sub(r'[\x00-\x1F]', '', repaired_str)  
        

        return json.loads(repaired_str)
    
    except json.JSONDecodeError as e:

        error_line = repaired_str.split('\n')[e.lineno - 1]
        error_ptr = ' ' * (e.colno - 1) + '^'
        print(f"parse failed in {e.lineno} line:\n{error_line}\n{error_ptr}")
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
                    raise ValueError("JSON parse failed")
            
def merge_fake_molecules(table):

    last_valid_key = None
    last_valid_entry = None
    
    for entry in table:
        
        for key in list(entry.keys()):
            if key == "XXXX":
                if last_valid_entry and last_valid_key:
                    
                    last_valid_entry[last_valid_key].extend(entry[key])
                del entry[key] 
            else:
                
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
    Calculate the charges and spin multiplicities of the molecular system
    Parameters:
        mol_input: Molecular SMILES expression or list of atomic symbols (e.g. ['C','H','H','H','H'])
        charge 
        is_atom
    Return:
        (charge, multiplicity)
    """

    if isinstance(atoms[0],int):
        total_e = sum(atoms) - charge
    else:
        total_e = sum(pt.GetAtomicNumber(symbol) for symbol in atoms) - charge
    multiplicity = 1  # default 

    if total_e % 2 == 1:  
        multiplicity = 2
    elif is_atom:  
        atomic_num = pt.GetAtomicNumber(atoms[0])
        valence_e = get_valence_electrons(atomic_num)
        multiplicity = 2 * (valence_e % 2) + 1 if valence_e < 14 else 3
    elif total_e == 0:  
        raise ValueError("The total number of electrons cannot be zero.")

    return (charge, multiplicity)

def generate_gauss_input_file(gjf_file,atoms,coords,charge,multiplicity,method,jobtype,cpu,memory,title="Generated by QuantumCalcAgent",ckpt=False):
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

    period = pt.GetPeriod(atomic_num)
    group = pt.GetGroup(atomic_num)
    
    if period == 1: return 1
    if period == 2: return min(group, 8)
    if period >=3:  
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

    element_counter = Counter()
    
    
    for sublist in nested_list:
        element_counter.update(sublist)
    
    return [element for element, count in element_counter.items() if count > threshold and len(element) > 5]

def parse_atom_coordinate(line):
    """
    Parse the atomic coordinate lines and extract the element symbols and 3D coordinates
    
    Parameters
    line (str): Atomic coordinate string, format example:
                "Rh-0.84075200 -1.17783800 -0.44030400"
                "Rh -0.84075200 -1.17783800 -0.44030400"
                "Rh-0.84075200 -1.17783800 -0."
    
    Return:
    tuple: (Element symbol, [x, y, z])
           or None (When the formats do not match)
    """

    pattern = re.compile(
        r'^\s*'                          
        r'([A-Za-z]{1,3})'               
        r'\d*'                           
        r'([-+]?\s*)?'                   
        r'([+-]?\d+\.?\d*)'              
        r'\s+([+-]?\d+\.?\d*)'           
        r'\s+([+-]?\d+\.?\d*)'           
        r'\s*$'                          
    )
    
    match = pattern.match(line)
    if not match:
        return None
    
    element = match.group(1)
    x_coord = match.group(3)
    if match.group(2):  
        x_coord = match.group(2).strip() + x_coord
    
    return element, [x_coord, match.group(4), match.group(5)]

def parse_atom_coordinatev2(line):
    """
    Parse the atomic coordinate lines and extract the element symbols and 3D coordinates
    
    Parameters
    line (str): Atomic coordinate string, format example: "N,0,0.5598740651,2.5536949817,1.2494624052"
    
    Return:
    tuple: (Element symbol, [x, y, z])
    """

    parts = line.strip().split(',')
    element = parts[0]
    xyz_list = [float(coord) for coord in parts[-3:]]
    
    return element, xyz_list