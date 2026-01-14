from ase.io import read, write
from ase.neighborlist import NeighborList
from ase.data import covalent_radii
import numpy as np


def grow_by_one_shell(is_selected, nl, safe_mask=None):
    """
    从当前选中原子出发，扩展一层邻近原子（不筛选元素）。
    is_selected: 布尔数组，当前已被选择的原子
    nl: ASE NeighborList（必须已 update）
    safe_mask: 布尔数组，安全元素位置为 True；这些原子即使在邻近壳层内也不选中
    """
    new_selected = is_selected.copy()
    seed_indices = np.where(is_selected)[0]

    for i in seed_indices:
        indices, offsets = nl.get_neighbors(i)
        for j in indices:
            if safe_mask is not None and safe_mask[j]:
                continue
            new_selected[j] = True

    if safe_mask is not None:
        new_selected[safe_mask] = False
    
    return new_selected


def select_shell_region(
        cif_file,
        seed_elems,
        safe_elems=None,
        n_shells=2,
        cutoff_scale=1.05,
        verbose=True):
    """
    参数:
        cif_file:  CIF 文件路径
        seed_elems: 初始选定元素集合，例如 {"Si","Al","P"}
        safe_elems: 安全元素集合；若原子属于安全元素，则不选中（包括种子与扩展壳层），例如 {"Na","K"}
        n_shells: 扩展壳层数 (int)
        cutoff_scale: covalent_radii 乘子的 cutoff，默认 1.1
        verbose: 是否打印过程信息

    返回:
        selected_atoms (ASE Atoms)
        remaining_atoms (ASE Atoms)
    """

    # 1. 读取 CIF
    atoms = read(cif_file)
    N = len(atoms)
    if verbose:
        print(f"[INFO] Loaded {cif_file}, atoms = {N}")

    safe_elems = set(safe_elems or [])
    safe_mask = np.array([atom.symbol in safe_elems for atom in atoms], dtype=bool)

    # 2. 初始种子区
    is_selected = np.array([atom.symbol in seed_elems for atom in atoms])
    if safe_elems:
        is_selected[safe_mask] = False
    if verbose:
        print(f"[INFO] Initial seed count = {is_selected.sum()}")

    # 3. 构建邻接
    cutoffs = [covalent_radii[a.number] * cutoff_scale for a in atoms]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    # 4. 扩展 n_shells
    for shell in range(1, n_shells + 1):
        is_selected = grow_by_one_shell(is_selected, nl, safe_mask=safe_mask if safe_elems else None)
        if verbose:
            print(f"[INFO] After shell {shell}, selected = {is_selected.sum()}")

    # 5. 输出区域
    selected_indices = np.where(is_selected)[0]
    remaining_indices = np.where(~is_selected)[0]

    selected_atoms = atoms[selected_indices]
    remaining_atoms = atoms[remaining_indices]

    if verbose:
        print(f"[INFO] Final selected = {len(selected_atoms)}, remaining = {len(remaining_atoms)}")

    return selected_atoms, remaining_atoms
