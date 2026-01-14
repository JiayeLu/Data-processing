import os
from pathlib import Path
from collections import defaultdict
from ase.io import read

# --- 用户定义区 ---
reactions = [
    'CH3OH-HMB', 'CH3OH-TMB', 'CH3OH-Toluene',
    'CH3Z-HMB', 'CH3Z-TMB', 'CH3Z-Toluene'
]
merged_dir = 'merged3'  # 存放合并后文件的文件夹
# ----------------

def check_atom_consistency():
    merged_path = Path(merged_dir)
    if not merged_path.exists():
        print(f"错误：目录 '{merged_dir}' 不存在。")
        return

    # stats 存储结构: {(zeo, react): {natoms: [filenames]}}
    stats = defaultdict(lambda: defaultdict(list))
    
    print(f"正在读取 {merged_dir} 中的文件并以同系列为基准校验原子数...")

    all_files = list(merged_path.glob("*.cif"))
    
    for f_path in all_files:
        f_name = f_path.name
        
        # 1. 匹配 reaction 列表
        found_react = next((r for r in reactions if r in f_name), None)
        if not found_react:
            continue

        try:
            # 2. 改进的 Zeo 提取逻辑：利用 -HZ- 关键字定位 
            # 示例: Scaling-ANO-O2-Mg-CH3OH-HMB-a_2.cif
            # 先去掉 Scaling- 前缀
            clean_name = f_name.replace("Scaling-", "")
            
            # 找到反应物部分的起始索引，截取其之前的内容
            react_idx = clean_name.find(found_react)
            # temp_str 此时为 "ANO-O2-Mg-"
            temp_str = clean_name[:react_idx].strip('-')
            
            # 将剩余部分按 '-' 拆分，最后一段是金属位点 (m_site)，前面是 zeo
            parts = temp_str.split('-')
            if len(parts) > 1:
                zeo = "-".join(parts[:-1])  # 结果: ANO-O2
                m_site = parts[-1]          # 结果: Mg
            else:
                zeo = parts[0]              # 针对没有连字符的单名
                m_site = "unknown"

            key = (zeo, found_react)
            
            # 3. 读取原子数 [cite: 51, 65]
            atoms = read(f_path)
            natoms = len(atoms)
            
            stats[key][natoms].append(f_name)
        except Exception as e:
            print(f"  [Error] 无法处理文件 {f_name}: {e}")

    # --- 报告输出 ---
    print("\n" + "="*65)
    print(f"{'拓扑结构':<15} | {'反应体系':<20} | {'状态'}")
    print("-" * 65)

    has_conflict = False
    sorted_keys = sorted(stats.keys())
    
    for key in sorted_keys:
        zeo, react = key
        count_dict = stats[key]
        
        if len(count_dict) == 1:
            n = list(count_dict.keys())[0]
            print(f"{zeo:<15} | {react:<20} | ✅ 一致 ({n} atoms)")
        else:
            has_conflict = True
            print(f"{zeo:<15} | {react:<20} | ❌ 冲突！")
            for n, files in count_dict.items():
                print(f"    -> 原子数 {n}: {len(files)} 个文件 (示例: {files[0]})")
    
    print("="*65)
    if not has_conflict:
        print("🎉 所有同类体系原子数完全一致！")
    else:
        print("⚠️ 发现原子数不一致，请核查合并时的去重或过滤逻辑。")

if __name__ == "__main__":
    check_atom_consistency()