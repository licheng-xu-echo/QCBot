import fitz,re,json,os,shutil
from rdkit import Chem
import numpy as np
from collections import Counter
from .valid import check_hydrogens,is_molecule_confortable
from copy import deepcopy
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

def extract_coords_info_from_pdf(pdf_files,dst_dir,cpu=32,mem=80):
    # match single column
    single_atom_pattern = re.compile(r'^\s*([A-Za-z]{1,3})\d*[-+]?\s*([+-]?\d+\.?\d*)\s+([+-]?\d+\.?\d*)\s+([+-]?\d+\.?\d*)\s*$') # fix Ru-0.11111 situation 
    # match double column
    dual_atom_pattern = re.compile(
        r'^\s*([A-Za-z]{1,3})\d*[-+]?\s*([+-]?\d+\.?\d*)\s+([+-]?\d+\.?\d*)\s+([+-]?\d+\.?\d*)\s*$'
        r'^\s*([A-Za-z]{1,3})\d*[-+]?\s*([+-]?\d+\.?\d*)\s+([+-]?\d+\.?\d*)\s+([+-]?\d+\.?\d*)\s*$')
    zc_pattern = re.compile(r"^[A-Z][a-z]?,-?\d+,-?\d+\.\d+,-?\d+\.\d+,-?\d+\.\d+$")

    for pdf_file in pdf_files:

        filename = os.path.basename(pdf_file)[:-4]

        #print(pdf_file)
        geom_dir = f"./{dst_dir}/{filename}"
        doc = fitz.open(pdf_file)
        text_lst = [[i+1,page.get_text()] for i,page in enumerate(doc)]


        # remove page index and empty lines
        text_wo_page_id = []
        for page_id,page_text in text_lst:
            page_lines = page_text.split("\n")
            #print(page_lines)
            page_lines = clean_empty_lines(drop_space(page_lines))
            text_wo_page_id.append(page_lines)
        header_and_footer_elements = find_header_and_footer_elements(text_wo_page_id,len(text_wo_page_id)-10)
        text_wo_page_id = [drop_page_with_specific_elements(page_lines,header_and_footer_elements) for page_lines in text_wo_page_id]
        text_wo_page_id = [drop_page_id_line(page_lines,page_id) for page_id,page_lines in enumerate(text_wo_page_id)]

        lines = []
        for page_idx,page_lines in enumerate(text_wo_page_id):
            lines.extend([(page_idx,line) for line in page_lines])
        match_text = []
        zc_match_text = []
        for idx,(page_idx,line) in enumerate(lines):
            single_match = single_atom_pattern.match(line)
            dual_match = dual_atom_pattern.match(line)
            zc_match = zc_pattern.match(line)
            if single_match or dual_match:
                match_text.append([idx,page_idx,line])
            if zc_match:
                zc_match_text.append([idx,page_idx,line])
        if len(match_text) == 0 and len(zc_match_text) == 0:
            print(f"[WARN] No match found in {filename}, skip")
            continue
        elif len(match_text) > 0:
            match_text_groups = split_into_continuous_groups(match_text)
            
            title_lines = []
            for group in match_text_groups:
                start_idx = group[0][0]-10
                end_idx = group[0][0]
                title_lines.append("//".join([line[1] for line in lines[start_idx:end_idx]]))
            
            
            geom_inf = [[parse_atom_coordinate(line[2]) for line in group] for group in match_text_groups]
        else:
            match_text = deepcopy(zc_match_text)
            match_text_groups = split_into_continuous_groups(match_text)
            
            title_lines = []
            for group in match_text_groups:
                start_idx = group[0][0]-10
                end_idx = group[0][0]
                title_lines.append("//".join([line[1] for line in lines[start_idx:end_idx]]))
            
            
            geom_inf = [[parse_atom_coordinatev2(line[2]) for line in group] for group in match_text_groups]

        if len(geom_inf) > 0:
            atoms = [[line[0] for line in group] for group in geom_inf]
            coords = [np.array([list(map(float,line[1])) for line in group]) for group in geom_inf]

            xyz_dir = f"{geom_dir}/xyz"
            os.makedirs(xyz_dir,exist_ok=True)
            gjf_dir = f"{geom_dir}/gjf"
            os.makedirs(gjf_dir,exist_ok=True)
            skip_whole_file = False
            for i,(atom,coord,title) in enumerate(zip(atoms,coords,title_lines)):
                # to handle the abnormal cases of structure extraction
                if len(atom) == 1:
                    print(f"[WARN] Only one atom found, skip. File: {filename} molecular index: {i}")
                    continue
                elif not "H" in atom:
                    print(f"[WARN] No hydrogen atom found, skip. File: {filename} molecular index: {i}")
                    continue
                elif not is_molecule_confortable(atom,coord,threshold_ratio=0.5):
                    print(f"[WARN] Molecule is not comfortable, skip. File: {filename} molecular index: {i}")
                    skip_whole_file = True
                    break
                elif not check_hydrogens(atom,coord,threshold_dist=2.0):
                    print(f"[WARN] There are no heavy atom around hydrogen atoms, skip. File: {filename} molecular index: {i}")
                    skip_whole_file = True
                    break

                xyz_file = f"{xyz_dir}/{i}.xyz"
                gjf_file = f"{gjf_dir}/{i}.gjf"
                charge,multi = get_charge_and_mult(atom,charge=0)
                symbol_pos_to_xyz_file(atom,coord,xyz_file,title=title)
                generate_gauss_input_file(gjf_file,atom,coord,charge,multi,method="p b3lyp/def2svp geom=PrintInputOrient",jobtype="opt",cpu=cpu,memory=mem,title=title)
            if skip_whole_file:
                print(f"[WARN] ++++++++++++++++++++ Skip whole file: {filename} ++++++++++++++++++++")
                shutil.rmtree(geom_dir)
            print(f"[INFO] ~~~~~~~~~~~~~~~~~~~~~~~~~~ {len(match_text_groups)} geometries found in {filename}, gjf and xyz files generated ~~~~~~~~~~~~~~~~~~~~~~~~~~")


def assign_file_type(xyz_files):
    filtered_words = ["Units","units",
    "Respresents","respresents",
    "Exhibits","exhibits",
    "Results","results",
    "Outside","outside",
    "Points","points",
    "Reactants","reactants",
    "Products","products",
    "Reagents","reagents",
    "Components","components"]
    int_filtered_words = ["POINT","QUINT","point","quint"]
    ts_geom = []
    int_geom = []
    other_geom = []
    for xyz_file in xyz_files:
        with open(xyz_file,"r") as fr:
            lines = fr.readlines()
        title = lines[1].strip()
        title_blks = title.split("//")
        match_flag = False
        for title_blk in title_blks[-3:]:
            if "TS" in title_blk or 'ts-' in title_blk or '-ts' == title_blk[-3:] or 'Ts-' in title_blk:
                ts_geom.append((xyz_file,title_blk))
                match_flag = True
                break
            elif "INT" in title_blk or "IM" in title_blk or 'IN' in title_blk or 'Int' in title_blk or 'int-' in title_blk or title_blk.isdigit():
                filtered_word_match_flag = False
                for filtered_word in int_filtered_words:
                    if filtered_word in title_blk:
                        filtered_word_match_flag = True
                        break
                if not filtered_word_match_flag:
                    int_geom.append((xyz_file,title_blk))
                    match_flag = True
                    break
            

        if not match_flag:
            other_geom.append((xyz_file,title_blk))

    return ts_geom,int_geom,other_geom
