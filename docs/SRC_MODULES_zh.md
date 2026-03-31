# src 模块说明（MacroMapFF）

本文档给出 `src/macromapff` 下每个 Python 模块的职责说明，便于快速定位代码入口与数据流。

## 包入口与编排

- `src/macromapff/__init__.py`
  - 包元信息与版本导出。
- `src/macromapff/cli.py`
  - 命令行入口，负责解析子命令并调用工作流接口。
- `src/macromapff/workflow.py`
  - 高层编排逻辑：发现样本、构建数据库、调用 parameterize 产出 LAMMPS 数据。

## pipeline 层

- `src/macromapff/pipeline/__init__.py`
  - pipeline 包声明与迁移说明。
- `src/macromapff/pipeline/env_build.py`
  - 从结构文件 + LAMMPS 数据构建原子环境特征，输出每模块 `*_atom_env.csv`。
- `src/macromapff/pipeline/keymap_build.py`
  - 合并多模块 atom env，生成全局 `final_env_keymap.csv` 与合并日志。
- `src/macromapff/pipeline/hop_build.py`
  - 从 final keymap 构建 hop2/hop1/hop0 回退映射表。
- `src/macromapff/pipeline/multi_extract.py`
  - 从单模块 LAMMPS 拓扑提取 bond/angle/dihedral/improper 观测项，输出 `*_multiatom_observed.csv`。
- `src/macromapff/pipeline/multi_build.py`
  - 合并多模块 multiatom 观测数据，映射到全局 key_type，输出 `multiatom_master_keytype.csv`。
- `src/macromapff/pipeline/lammps_gen.py`
  - parameterize 阶段主程序：读数据库并为目标分子生成参数化 LAMMPS data 文件。

## pipeline/core 公共能力

- `src/macromapff/pipeline/core/__init__.py`
  - core 子包声明。
- `src/macromapff/pipeline/core/env.py`
  - 环境键（env_key）规范化、拆分与结构化索引相关工具。
- `src/macromapff/pipeline/core/atom_match.py`
  - 原子级匹配逻辑：按 env_key 精确匹配并执行 hop 回退。
- `src/macromapff/pipeline/core/multi_match.py`
  - 多体项匹配逻辑：按 interaction kind + key_type_tuple 匹配参数。
- `src/macromapff/pipeline/core/lammps_parse.py`
  - LAMMPS 数据文件解析能力（Masses、Atoms、拓扑等）。
- `src/macromapff/pipeline/core/lammps_write.py`
  - 将匹配结果写回 LAMMPS data 文件格式。

## 典型调用链

- build-db:
  - `cli.py` -> `workflow.py` -> `env_build.py` -> `multi_extract.py` -> `keymap_build.py` -> `hop_build.py` -> `multi_build.py`
- parameterize:
  - `cli.py` -> `workflow.py` -> `lammps_gen.py` -> `core/atom_match.py` + `core/multi_match.py` + `core/lammps_write.py`
