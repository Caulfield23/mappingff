# src 模块说明（MacroMapFF）

本文档给出 `src/macromapff` 下每个 Python 模块的职责说明，便于快速定位代码入口与数据流。

## 包入口与编排

- `src/macromapff/__init__.py`
  - 包元信息与版本导出。
- `src/macromapff/cli.py`
  - 单文件入口：包含 CLI 解析 + `Workflow` 编排类。
  - `Workflow.discover_samples(...)` 作为类内方法负责样本发现与输入资产路径匹配。
  - 多 sample 场景（build-db / add-samples）的编排职责全部在 `cli.py`。

## domain 层

- `src/macromapff/domain/__init__.py`
  - domain 统一导出入口，外层优先从该入口导入核心逻辑能力。

- `src/macromapff/domain/env.py`
  - env_key 规范化与拆分工具。
- `src/macromapff/domain/env_features.py`
  - 基于 RDKit 的原子环境特征提取。
- `src/macromapff/domain/atom_typing_core.py`
  - 原子类型匹配的纯逻辑。
- `src/macromapff/domain/term_enumeration.py`
  - 枚举 bonds/angles/dihedrals/impropers。
- `src/macromapff/domain/multiatom_match_core.py`
  - 多体匹配纯逻辑与 type 映射构建。
- `src/macromapff/domain/multiatom_observed.py`
  - 多体观测映射构建逻辑。
- `src/macromapff/domain/keymap_merge.py`
  - env_key 合并与统计逻辑。
- `src/macromapff/domain/multiatom_master_merge.py`
  - 多体主库合并逻辑。

## io 层

- `src/macromapff/io/input.py`
  - 仅负责输入读取与解析：结构读取、LAMMPS 解析、数据库 CSV 加载。
- `src/macromapff/io/output.py`
  - 仅负责输出写出：atom_env、keymap、multiatom、LAMMPS data 等结果文件。
- `src/macromapff/io/log.py`
  - 仅负责日志输出：build/missing、keymap merge、多体匹配与冲突日志。

## pipeline 层

- `src/macromapff/pipeline/__init__.py`
  - pipeline 统一导出入口，workflow 通过该入口组织流程组件。
- `src/macromapff/pipeline/atommap_sample.py`
  - 只负责单 sample 的 `*_AtomMap.csv` 构建。
- `src/macromapff/pipeline/keymap_hop.py`
  - 合并多模块 env_key，输出 `Global_AtomMap.csv`，并构建 hop2/hop1/hop0 回退映射表。
- `src/macromapff/pipeline/bondedterms_sample.py`
  - 只负责单 sample 的 `*_BondedTerms.csv` 提取。
- `src/macromapff/pipeline/multiatom_master.py`
  - 合并多模块多体观测，输出 `Global_BondedTerms.csv`。
- `src/macromapff/pipeline/parameterize.py`
  - parameterize 阶段主程序，生成参数化 LAMMPS data（含原子匹配与多体匹配流程封装）。

## 典型调用链

- 统一导入约定：
  - IO 能力统一从 `macromapff.io` 导入。
  - 编排能力统一从 `macromapff.pipeline` 导入。
  - 纯逻辑能力统一从 `macromapff.domain` 导入（优先）。

- build-db:
  - `cli.py(Workflow)` -> `atommap_sample.py` -> `bondedterms_sample.py` -> `keymap_hop.py` -> `multiatom_master.py`
- add-samples:
  - `cli.py(Workflow)` 仅对新样本运行 `atommap_sample.py` + `bondedterms_sample.py`
  - 输出目录直接为 `db-dir/*_env`，与已有 sample env 同级
  - 若 module 同名则覆盖对应 `*_env`
  - 通过扫描 `db-dir` 下全部 `*_env` 自动重建 `Global_AtomMap.csv` / `hop_env/*` / `Global_BondedTerms.csv`
- parameterize:
  - `cli.py(Workflow)` -> `parameterize.py` + `io/output.py`
