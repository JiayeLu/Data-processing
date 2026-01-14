import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import NeighborList, neighbor_list
from ase.data import covalent_radii
from pathlib import Path

def extract_molecule_from_guest(guest_atoms, seed_elems={"Si", "Al", "P", "Ge", "Mg", "Zn"}):
    """
    从带有骨架的 Guest 文件中提取反应物分子
    逻辑：识别骨架原子并剔除，但绝对保护 C 原子
    """
    atoms = guest_atoms.copy()
    N = len(atoms)
    
    # 1. 构建邻居列表用于壳层扩展
    cutoffs = [covalent_radii[a.number] * 1.05 for a in atoms]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    # 2. 标记骨架种子 (Si, Al, P 等)
    is_skeleton_seed = np.array([atom.symbol in seed_elems for atom in atoms])
    
    # 3. 标记必须保护的原子 (所有 C 原子)
    is_protected = np.array([atom.symbol == "C" for atom in atoms])
    
    # 4. 扩展骨架区域：寻找与 Si/Al/P 相连的 O (通常扩展 1-2 层即可覆盖骨架氧)
    # 第一层扩展：Si -> O
    is_skeleton = is_skeleton_seed.copy()
    for i in np.where(is_skeleton_seed)[0]:
        indices, _ = nl.get_neighbors(i)
        for j in indices:
            # 如果邻居不是 C，则标记为骨架
            if not is_protected[j]:
                is_skeleton[j] = True
                
    # 5. 最终检查：从“骨架候选”中剔除 C 原子及其紧密相连的原子（双重保险）
    # 这样剩下的就是真正的反应物分子（C, H, 以及不属于骨架的 O/N 等）
    keep_mask = ~is_skeleton
    # 强制保留所有 C 及其连接的原子
    for i in np.where(is_protected)[0]:
        keep_mask[i] = True
        indices, _ = nl.get_neighbors(i)
        for j in indices:
            keep_mask[j] = True

    return atoms[keep_mask]

def merge_to_host(host, guest_molecule, tol_same=1.20, min_dist=0.30):
    """
    使用之前优化的逻辑，将提取出的分子安全合并到 Host 中
    """
    host_c = host.copy()
    guest_c = guest_molecule.copy()
    
    # 统一晶格
    guest_c.set_cell(host_c.cell, scale_atoms=False)
    guest_c.set_pbc(host_c.pbc)
    guest_c.wrap()

    # 去重逻辑 (元素相同且距离近则删除)
    combined = host_c + guest_c
    symbols = np.array(combined.get_chemical_symbols())
    i, j, d = neighbor_list('ijd', combined, cutoff=tol_same, self_interaction=False)
    
    mask_host = i < len(host_c)
    mask_guest = j >= len(host_c)
    mask_same_symbol = (symbols[i] == symbols[j])
    
    drop_indices = j[mask_host & mask_guest & mask_same_symbol]
    guest_indices_to_drop = np.unique(drop_indices) - len(host_c)
    
    keep_mask = np.ones(len(guest_c), dtype=bool)
    if len(guest_indices_to_drop) > 0:
        keep_mask[guest_indices_to_drop] = False
    
    final_atoms = host_c + guest_c[keep_mask]
    return final_atoms

def process_all_files(host_dir, guest_dir, output_dir, reactions):
    host_paths = list(Path(host_dir).glob("*.cif"))
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for h_path in host_paths:
        f_name = h_path.name  # Scaling-ANO-O2-Mg-HZ-a_2.cif
        
        # 1. 提取 Reaction 和 Zeo 名
        found_react = next((r for r in reactions if r in f_name), None)
        if not found_react:
            # 如果 host 里不含反应名，则根据你的逻辑匹配对应的 Guest
            # 假设你对每个反应都要跑一遍这个 host
            current_reactions = reactions 
        else:
            current_reactions = [found_react]

        for r in current_reactions:
            # 解析 Zeo 名：提取 Scaling- 和 -{Metal} 之间的所有部分
            # 针对 Scaling-ANO-O2-Mg-HZ-a_2.cif -> ANO-O2
            parts = f_name.replace("Scaling-", "").split("-")
            # 此时 parts 可能为 ['ANO', 'O2', 'Mg', 'HZ', 'a_2.cif']
            # 金属位点通常在 HZ 之前
            try:
                hz_idx = parts.index("HZ")
                m_site = parts[hz_idx - 1]
                zeo = "-".join(parts[:hz_idx - 1])
            except ValueError:
                continue

            guest_name = f"Scaling-{zeo}-Si-{r}-a_2.cif"
            g_path = Path(guest_dir) / guest_name
            
            if not g_path.exists():
                continue

            # 2. 执行处理
            host_atoms = read(h_path)
            # 去除 Host 中的 H
            non_h_indices = [atom.index for atom in host_atoms if atom.symbol != 'H']
            host_atoms = host_atoms[non_h_indices]
            
            guest_raw = read(g_path)
            # 第一步：只从 Guest 中提取反应物分子
            molecule = extract_molecule_from_guest(guest_raw)
            # 第二步：合并到 Host
            final_structure = merge_to_host(host_atoms, molecule)
            
            # 3. 保存
            out_name = f"Scaling-{zeo}-{m_site}-{r}-a_2.cif"
            write(output_dir / out_name, final_structure)
            print(f"Successfully processed: {out_name}")
if __name__ == "__main__":
    host_directory = "data"
    guest_directory = "Scaling-Si-6R"
    output_directory = "merged3"
    reactions_list = ["CH3Z-HMB", "CH3OH-TMB", "CH3Z-Toluene",
                      "CH3OH-HMB", "CH3OH-Toluene", "CH3Z-TMB"]

    process_all_files(host_directory, guest_directory, output_directory, reactions_list)
