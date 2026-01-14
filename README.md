ZeoBuilder: 自动化金属取代分子筛建模与校验工具ZeoBuilder 是一套专为计算化学设计的 Python 脚本工具集。它能够高效处理金属取代分子筛（Zeolite）模型的自动化搭建、区域选择及原子数校验，确保不同取代环境下的模型严格对齐。🌟 核心特性智能分子提取：从带有骨架的 Guest 模板中提取反应物，通过邻近壳层算法精准识别并剔除旧骨架。元素保护机制：在合并时通过元素符号匹配保护碳原子（C）及其官能团，防止误删。复杂拓扑识别：完美解析包含连字符的拓扑结构名（如 ANO-O2），支持动态锚点定位。区域选择扩展：支持基于原子种子的壳层扩展逻辑，方便提取局部活性位点结构。严格基准校验：强制以 Si-6R 体系作为“黄金标准”进行原子数校验。🛠️ 脚本功能清单1. merge_atoms_ultra_safe.py (模型构建)核心建模工具。读取 Host（金属取代骨架）和 Guest（含反应物的 Si 体系），在 $1.20$ Å 范围内识别并剔除重复骨架，同时保留反应物分子。2. select_frame.py (区域选择)利用 NeighborList 算法，从指定的种子原子（如 Si, Al, P）出发，向外扩展 $N$ 层壳层（Shells）。常用于从庞大的分子筛骨架中切取局部簇模型（Cluster Model）进行精细计算。3. check_atom_consistency.py (质量检查)自动化审计工具。扫描 merged 文件夹，解析 Zeolite 类型和 Reaction 类型，并与 Scaling-Si-6R 目录下的原始文件对比原子数，确保计算基准完全一致。🚀 快速开始1. 安装依赖Bashpip install ase numpy
2. 数据准备确保目录结构如下：test2/: 存放 Host 骨架code-test/: 存放 Guest 模板Scaling-Si-6R/: 原子数标准参考库merged/: 输出目录3. 运行工作流Python# 1. 执行合并
python merge_atoms_ultra_safe.py

# 2. (可选) 执行区域切取
python select_frame.py

# 3. 校验原子数
python check_atom_consistency.py
📂 文件命名规范项目采用以下动态解析逻辑，支持复杂命名：Scaling-{Zeolite}-{Metal}-{Reaction}-a_2.cifZeolite: 拓扑名，支持 AEI, ERI, ANO-O2 等 1。Metal: 取代金属，如 Mg, Zn, Al, Ga, Si 2。Reaction: 反应体系，如 CH3OH-HMB, CH3Z-Toluene 3。
