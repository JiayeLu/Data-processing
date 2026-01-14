# ZeoBuilder: 自动化金属取代分子筛建模与校验工具



**ZeoBuilder** 是一套专为计算化学设计的 Python 脚本工具集，旨在解决分子筛（Zeolite）模型构建中繁琐的重复劳动。它能够高效处理金属取代模型的自动化搭建、局部区域提取及原子数一致性校验。

---

## 🌟 核心特性

- [cite_start]**智能分子提取**：自动从带有骨架的 Guest 模板中提取反应物分子，通过邻近壳层算法精准剔除旧骨架 [cite: 53, 56, 60]。
- [cite_start]**元素级保护机制**：在合并过程中，通过元素符号匹配强制保护碳原子（C）及其连接的反应物官能团，防止因距离判定导致的误删 [cite: 55, 56]。
- **复杂拓扑识别**：采用动态锚点定位逻辑，完美解析包含连字符的复杂拓扑结构名（如 `ANO-O2`）。
- **灵活区域提取**：支持基于原子种子的壳层扩展算法（NeighborList），方便从庞大骨架中切取局部簇模型（Cluster Model）。
- **严格基准校验**：强制以 `Si-6R` 体系作为“黄金标准”进行原子数对齐校验，确保 Scaling Relations 分析的能量基准一致性。

---

## 🛠️ 脚本功能清单

### 1. `merge_atoms_ultra_safe.py` (模型构建)
[cite_start]核心建模脚本。它读取 Host（金属取代骨架）和 Guest（含反应物的 Si 体系模板），在 $1.20$ Å 的范围内识别并剔除重复骨架原子 [cite: 1, 35, 56]。
- [cite_start]**功能点**：支持 MIC (Minimum Image Convention) 周期性边界处理 [cite: 51, 52]。

### 2. `select_frame.py` (区域选择)
[cite_start]局部区域切取工具。利用 `NeighborList` 算法，从指定的种子原子（如 Si, Al, P）出发，向外扩展 $N$ 层壳层（Shells） [cite: 35, 42, 50]。常用于制备 QM/MM 或者是局部簇模型。

### 3. `check_atom_consistency.py` (质量检查)
[cite_start]自动化审计工具。扫描 `merged` 文件夹，解析 Zeolite 类型和 Reaction 类型，并与 `Scaling-Si-6R` 目录下的原始文件对比原子数 。

---

## 🚀 快速开始

### 1. 安装依赖
```bash
pip install ase numpy
