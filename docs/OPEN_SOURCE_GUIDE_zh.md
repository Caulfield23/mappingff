# MacroMapFF 开源化与工程管理文档

## 1. 文档目标

本文档定义 MacroMapFF 的开源仓库组织方式、工作流、版本策略、数据资产管理规则和迁移说明。

目标是让项目从“可跑”升级到“可维护、可协作、可发布”。

## 2. 项目定位

MacroMapFF 采用“Python 包 + CLI”双入口模式：

- Python 包：适合研究复用、函数级集成、Notebook 调用、后续二次开发。
- CLI：适合生产批处理、可复现实验、脚本化运行和 CI 自动化。

这比仅做脚本集合更利于长期维护，也比仅做 API 更方便非开发用户。

## 3. 命名和版本

- 分发名：MacroMapFF
- import 名：macromapff
- 命令名：MacroMapFF
- 当前版本：0.1.0

语义化版本建议：

- MAJOR：破坏性变更（参数格式/CLI 行为变化）
- MINOR：新增功能且兼容
- PATCH：缺陷修复，不改变行为契约

## 4. 目录结构说明

```
MacroMapFF/
  .github/workflows/
    ci.yml
    release.yml
  docs/
    OPEN_SOURCE_GUIDE_zh.md
  examples/
    ps_odms7poss_legacy/
      PS-oDMS7POSS.mol
      PS-oDMS7POSS.pdb
      segment1/
      segment2/
      segment3/
      segment4/
      outputs/
      scripts/                # legacy workflow, isolated from root source
  src/
    macromapff/
      __init__.py
      cli.py
      domain/
        __init__.py
        env.py
        env_features.py
        atom_typing_core.py
        term_enumeration.py
        multiatom_match_core.py
        multiatom_observed.py
        keymap_merge.py
        multiatom_master_merge.py
      io/
        __init__.py
        input.py
        output.py
        log.py
      pipeline/
        atom_env.py
        keymap_hop.py
        multiatom_observed.py
        multiatom_master.py
        parameterize.py
  tests/
  pyproject.toml
  README.md
```

说明：

- 核心算法逻辑放在 `src/macromapff/domain/`，仅处理内存对象，不做文件 IO。
- 文件解析/序列化放在 `src/macromapff/io/`，负责结构文件与 LAMMPS data 的读写。
- 流程编排集中在 `src/macromapff/pipeline/`，调用 domain + io 完成完整流程。
- 历史可运行工程（分子、segments、outputs、legacy scripts）已整体迁入 `examples/ps_odms7poss_legacy/`。
- 根目录仅保留开源代码与工程化配置，和 example 数据完全隔离。

## 5. 迁移策略（本次执行）

本次迁移遵循“只改工程组织，不改算法逻辑”：

1. 将历史 Python 脚本从 `scripts/` 迁移到 `src/macromapff/`，并按 domain/io/pipeline 分层归位。
2. 在 `examples/ps_odms7poss_legacy/scripts/` 下保留 legacy 启动器，转发到新包模块。
3. 新增 `pyproject.toml`，支持 `pip install -e .` 和 `MacroMapFF` 命令。
4. 新增 CI、Release、pre-commit、LICENSE、README。
5. 将历史输入/输出和 segment 资产打包到单一 example，实现与现有源码彻底隔离。

不变项：

- 核心参数匹配逻辑
- 数据库构建逻辑
- impropers/bonds/angles/dihedrals 的既有实现行为

## 6. 工作流规范

### 6.1 分支管理

- `main`：稳定分支
- `feature/*`：功能开发
- `fix/*`：缺陷修复
- `docs/*`：文档更新

### 6.2 提交流程

建议模板：

- feat: add xxx
- fix: correct xxx
- docs: update xxx
- refactor: move xxx without behavior change

### 6.3 Pull Request 规范

每个 PR 应包含：

- 变更目的
- 影响范围（逻辑/接口/文档）
- 验证方法（命令、输入、输出）
- 回滚方式（如适用）

## 7. CI / 发布策略

### 7.1 CI（ci.yml）

触发：push main、pull request

步骤：

1. 安装包（editable）
2. CLI smoke test（`MacroMapFF --help`）
3. 运行 `pytest`

### 7.2 Release（release.yml）

触发：`v*` tag

产物：wheel + sdist（作为可下载构建产物）

后续可扩展到 PyPI 发布和 Release Note 自动生成。

## 8. 数据与示例管理

### 8.1 当前策略

- canonical legacy example 统一放在 `examples/ps_odms7poss_legacy/`
- 根目录不再放置历史输入输出资产
- `examples/ps_odms7poss_legacy/outputs/` 为历史生成数据快照

### 8.2 原则

- 输入样例尽量最小化且可复现
- 输出大文件不纳入 Git（建议用 artifact 或外部存储）
- 对应文档固定命令与预期结果

## 9. 运行模式

### 9.1 包命令模式（推荐）

仅保留顶层 CLI：

```bash
MacroMapFF build-db --help
MacroMapFF add-samples --help
MacroMapFF parameterize --help
```

### 9.2 兼容脚本模式（保持可用）

```bash
cd examples/ps_odms7poss_legacy/scripts
./rebuild_databases.sh
./generate_from_current_db.sh
```

## 10. 后续里程碑建议

- v0.1.1：补全测试样例与黄金输出比对
- v0.2.0：对 example 增加最小化子集与自动下载机制（可选）
- v0.3.0：稳定 API 层（`macromapff.api`）
- v1.0.0：冻结 CLI 接口契约、发布可引用论文/文档

## 11. 维护者检查清单

每次发布前：

1. `pip install -e .` 成功
2. `MacroMapFF --help` 成功
3. `scripts/rebuild_databases.sh` 可运行
4. `scripts/generate_from_current_db.sh` 可运行
5. `examples/ps_odms7poss_legacy/scripts/rebuild_databases.sh` 可运行
6. `examples/ps_odms7poss_legacy/scripts/generate_from_current_db.sh` 可运行
7. 示例输入存在，文档链接有效
8. CI 全绿

## 12. 兼容性声明

本次仅做开源工程化重构，不修改算法逻辑；若观察到结果变化，应优先检查：

- 运行路径和相对路径
- Python 环境依赖版本
- 输入输出文件是否一致

## 13. generate 映射数据流（详细）

本节描述当前实现中，从数据库构建到新分子 parameterize 的完整数据流，重点解释 atom 与 multiatom 的匹配过程，以及跨表查询如何工作。

### 13.1 总览（两阶段）

阶段 A：build-db（构建数据库）

1. 样本发现：扫描 samples 目录下所有 .lmp，并配对结构文件（.mol/.mol2/.pdb）。
2. 单模块原子环境提取：为每个模块输出 {module}_AtomMap.csv。
3. 单模块多体观测提取：输出 {module}_BondedTerms.csv。
4. 全模块环境合并：输出 Global_AtomMap.csv（包含 global_type_ids 多值映射列）。
5. hop 回退库构建：输出 hop_env/hop2_KeyMap.csv、hop_env/hop1_KeyMap.csv、hop_env/hop0_KeyMap.csv。
6. 多体主库构建：输出 Global_BondedTerms.csv。

阶段 B：parameterize（generate 新分子映射并写出 LAMMPS）

1. 载入 hop2/hop1/hop0 与 multiatom_master。
2. 对新分子每个原子做环境匹配，得到 atom 的全局 key_type。
3. 枚举 bonds/angles/dihedrals/impropers，做多体匹配。
4. 将匹配结果写入最终 LAMMPS data 文件。

### 13.2 build-db 的数据流细节

#### 13.2.1 {module}_AtomMap.csv（atom_env）

- 输入：模块结构文件 + 对应 .lmp。
- 对每个 atom 生成 env_key（当前最多 hop2），并保留该 atom 在样本中的 OPLS type 与 LJ 参数。
- 输出核心字段：
  - module, atom_index, atom_name
  - opls_type_id, opls_type_name
  - charge, sigma, epsilon
  - env_key

该表是后续 keymap 合并与 multiatom_observed 的共同输入。

#### 13.2.2 Global_AtomMap.csv（keymap_hop）

- 读取所有模块 atom_env 行，并按 canonical env_key 合并。
- 生成全局 key_id（稳定排序后编号）。
- Global_AtomMap.csv 保存 key 级别环境与均值参数。
- 其中第二列 ``global_type_ids`` 保存 type 到 key_id 的桥接关系（多值，分号分隔），
  每个值格式为 ``module_xx``。

这张 type_stats 表是后续 multiatom_master 的关键跨表桥。

#### 13.2.3 hop2/hop1/hop0_KeyMap.csv（keymap_hop）

- 从 Global_AtomMap.csv 聚合得到三个粒度的回退库。
- 每行包含：
  - source_key_ids（该聚合行由哪些 key_id 合并而来）
  - charge_mean, sigma_mean, epsilon_mean, mass_mean
  - env_key
  - 环境拆分列（hop1_shell, hop2_shell 位于末尾）

其中 hop0 的 source_key_ids 在 multiatom_master 中用于构建 key 等价类（见 13.2.4）。

#### 13.2.4 Global_BondedTerms.csv（multiatom_master）

该步骤发生两次跨表映射：

第一次跨表：模块 type -> 全局 key_id

- 来源表：Global_AtomMap.csv 的 ``global_type_ids`` 列
- 查询键：``segmentX_xx``（由模块名和 opls_type_id 拼接）
- 结果：key_id

第二次跨表：key_id -> key 等价类（key_type slot）

- 来源表：hop0_KeyMap.csv
- 用 source_key_ids 通过并查集（DSU）求连通分量
- 结果：每个 key_id 映射到一个 key class（如 [12, 29, 31]）

最后将模块观测表中的 lmp_type_tuple 转成 key_type_tuple（每个位置是允许 key 集合），并按 interaction_kind + key_type_tuple 合并，得到 Global_BondedTerms.csv。

### 13.3 parameterize(generate) 的数据流细节

#### 13.3.1 原子项匹配（atom_typing_core）

输入库：hop2_KeyMap.csv、hop1_KeyMap.csv、hop0_KeyMap.csv。

预处理：每个 hop 库都建两套索引

1. env 精确索引
  - 键：env_key（整串 JSON）
  - 值：参数与 key_ids

2. structured 索引
  - 键：拆分列 tuple（z/formal_charge/.../hop1_shell/hop2_shell/neighbor_sig/bond_kinds）
  - 值：参数与 key_ids

对新分子每个 atom 的查找顺序：

1. 计算该 atom 的 env 特征（hop 深度 2）
2. 先查 hop2
  - 先 env 精确命中
  - 再 structured 命中
3. 若未命中，按 fallback 顺序查 hop1 -> hop0
4. 命中后得到：
  - global_key_id（默认取 key_ids 第一个）
  - global_key_ids（该 atom 的候选 key 集）
  - charge/sigma/epsilon/mass

输出：

- atom_records（供写 Atoms 与后续多体匹配）
- atom_index_key_types.csv（人工审查用）

#### 13.3.2 多体项匹配（multiatom_match_core）

输入库：Global_BondedTerms.csv。

读取后构建倒排索引：

- 可逆项（bond/angle/dihedral）：按 (位置, key_type) 建倒排
- improper：按中心位 key_type 建倒排

匹配时步骤：

1. 先从结构枚举 terms：
  - bonds
  - angles
  - dihedrals
  - impropers
2. 对每个 term，读取每个原子的 global_key_ids，形成 key 选项集合。
3. 对选项做笛卡尔积，得到候选 key_tuple。
4. 用倒排索引筛候选 pattern，再做 slot 级集合包含判断：
  - 可逆项同时检查正向与反向
  - improper 检查中心位固定 + 其余三位排列匹配
5. 若命中多个候选：
  - 记录 ambiguous 日志
  - 采用稳定排序后的第一组 key_tuple 与第一组 coeff 作为最终选择
6. 若未命中：
  - 在当前配置下记 WARN，不中断（strict_missing=False）

#### 13.3.3 写出 LAMMPS（io/output）

- atom 的局部 type 编号：由已使用的 global_key_id 去重后重映射为 1..N。
- multiatom 的 type 编号：按 coeff 去重映射。
- 最终写出完整 LAMMPS data：Masses、Pair Coeffs、Atoms、Bonds、Angles、Dihedrals、Impropers 及对应 Coeffs。

### 13.4 跨表查询键一览（速查）

1. 模块原子类型 -> 全局 key_id
  - 表：Global_AtomMap.csv
  - 键：global_type_ids 中的 ``segmentX_xx``

2. 全局 key_id -> 等价 key class
  - 表：hop0_KeyMap.csv
  - 键：source_key_ids 连通关系

3. 新分子 atom 环境 -> 原子参数/候选 key_ids
  - 表：hop2/hop1/hop0_KeyMap.csv
  - 键：先 env_key，后 structured tuple

4. term 的 key_tuple -> 多体 coeff
  - 表：Global_BondedTerms.csv
  - 键：interaction_kind + key_type_tuple（slot 集合匹配）

### 13.5 人工核查推荐顺序

1. 先看 parameterize.log：missing/ambiguous/cache 命中率。
2. 看 atom_index_key_types.csv：候选 key 是否合理。
3. 看 Global_AtomMap.csv 与 hop2/hop1/hop0：回退层是否符合预期。
4. 看 Global_BondedTerms.csv：关键 interaction 的 key_type_tuple 与 coeff 是否覆盖。
