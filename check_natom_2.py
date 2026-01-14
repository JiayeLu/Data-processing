import os
from pathlib import Path
from ase.io import read

# --- 用户定义区 ---
reactions = [
    'CH3OH-HMB', 'CH3OH-TMB', 'CH3OH-Toluene',
    'CH3Z-HMB', 'CH3Z-TMB', 'CH3Z-Toluene'
]
merged_dir = Path('merged3')           # 待检查的合并后目录
ref_dir = Path('Scaling-Si-6R')       # 绝对参考标准目录 (Si 体系)
# ----------------

def check_against_si_reference():
    if not merged_dir.exists() or not ref_dir.exists():
        print(f"错误：请检查路径是否存在：\nMerged: {merged_dir}\nReference: {ref_dir}")
        return

    print(f"开始校验：所有结构必须与 {ref_dir} 中的 Si 模板原子数一致...\n")
    
    results = {
        "pass": 0,
        "fail": [],
        "missing_ref": []
    }

    # 遍历合并后的文件
    for f_path in merged_dir.glob("*.cif"):
        f_name = f_path.name
        
        # 1. 识别 reaction 锚点
        found_react = next((r for r in reactions if r in f_name), None)
        if not found_react:
            continue

        try:
            # 2. 动态提取 zeo 拓扑名 (支持 ANO-O2 等含连字符名称)
            # 逻辑：去除 Scaling- 前缀，截取到 reaction 之前，再去掉最后的金属位点
            clean_name = f_name.replace("Scaling-", "")
            react_idx = clean_name.find(found_react)
            temp_str = clean_name[:react_idx].strip('-')
            
            zeo_parts = temp_str.split('-')
            if len(zeo_parts) > 1:
                zeo = "-".join(zeo_parts[:-1]) # 结果示例: ANO-O2
            else:
                zeo = zeo_parts[0]              # 结果示例: AEI

            # 3. 构建对应的 Si 参考文件名
            # 格式: Scaling-{zeo}-Si-{reaction}-a_2.cif
            ref_name = f"Scaling-{zeo}-Si-{found_react}-a_2.cif"
            ref_path = ref_dir / ref_name

            if not ref_path.exists():
                results["missing_ref"].append(f"{f_name} (缺少标准: {ref_name})")
                continue

            # 4. 读取并比对原子数
            atoms_merged = read(f_path)
            atoms_ref = read(ref_path)
            
            n_merged = len(atoms_merged)
            n_ref = len(atoms_ref)

            if n_merged == n_ref:
                results["pass"] += 1
            else:
                results["fail"].append({
                    "file": f_name,
                    "n": n_merged,
                    "ref_n": n_ref,
                    "ref_file": ref_name
                })

        except Exception as e:
            print(f"  [Error] 处理 {f_name} 时发生异常: {e}")

    # --- 输出最终报告 ---
    print("="*70)
    print(f"{'拓扑-反应体系':<40} | {'校验结果'}")
    print("-" * 70)

    if results["fail"]:
        for item in results["fail"]:
            print(f"❌ {item['file']:<38} | 当前:{item['n']} vs 标准:{item['ref_n']}")
    
    if results["missing_ref"]:
        print("-" * 70)
        print("以下文件因缺少 Si 参考模板被跳过:")
        for m in results["missing_ref"]:
            print(f"  - {m}")

    print("="*70)
    print(f"📊 统计汇总:")
    print(f"   ✅ 通过数量: {results['pass']}")
    print(f"   ❌ 失败数量: {len(results['fail'])}")
    print(f"   ⚠️ 缺失标准: {len(results['missing_ref'])}")
    print("="*70)

    if not results["fail"] and not results["missing_ref"]:
        print("🎉 完美！所有合并后的结构原子数均与 Si 参考标准一致。")

if __name__ == "__main__":
    check_against_si_reference()