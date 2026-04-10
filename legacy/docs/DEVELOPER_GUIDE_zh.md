# MacroMapFF 开发者指南

本文档面向开发者和维护者，说明内部架构、边界约束与模块职责。

## 1. 内部规则

- 函数优先：除非需要持续可变状态，否则优先使用函数。
- 边界清晰：
  - `cli`：只做参数解析与命令分发。
  - `pipeline`：只做流程编排与调用顺序组织。
  - `domain`：只做业务规则与匹配逻辑。
  - `io`：只做输入解析、输出写出与日志。
- 不做伪配置化：
  - 匹配策略、fallback 顺序、内部常量保持在代码内部。
  - 仅将用户运行时输入暴露为 CLI 参数。
- 不做无意义透传：
  - 参数只传给真正消费该参数的函数。

## 2. 包结构

- `src/macromapff/cli.py`：CLI 入口。
- `src/macromapff/pipeline/`：工作流编排层。
- `src/macromapff/domain/`：匹配与合并逻辑层。
- `src/macromapff/io/`：解析、写出、日志层。

## 3. 工作流架构

### `build-db`
1. 发现样本输入（结构文件与 `.lmp`）。
2. 生成样本原子映射（`*_AtomMap.csv`）。
3. 生成样本多体观测（`*_BondedTerms.csv`）。
4. 合并全局原子映射（`Global_AtomMap.csv`）。
5. 构建 hop 回退数据库（`hop2/hop1/hop0_KeyMap.csv`）。
6. 合并全局多体数据库（`Global_BondedTerms.csv`）。

### `add-samples`
- 追加/覆盖样本后，执行与 `build-db` 相同的全局重建流程。

### `parameterize`
1. 读取 hop 数据库与全局多体数据库。
2. 执行原子匹配（`hop2 -> hop1 -> hop0` 回退）。
3. 枚举并匹配多体项。
4. 写出最终 LAMMPS 文件与日志。

## 4. 模块职责

### 4.1 CLI
- `macromapff.cli`
  - `main`：命令行解析与分发。

### 4.2 Pipeline
- `pipeline.workflow`
  - `discover_samples`：发现样本输入。
  - `build_sample_envs`：构建样本 `_env` 产物。
  - `discover_sample_env_records`：收集合并记录。
  - `merge_database`：重建全局数据库产物。
  - `build_db`、`add_samples`、`parameterize`：顶层工作流。
- `pipeline.atommap_sample`
  - `build_sample_atommap`：生成单样本 atom map。
- `pipeline.bondedterms_sample`
  - `extract_bondedterms_sample`：生成单样本 bonded 观测表。
- `pipeline.keymap_hop`
  - `build_keymap`：合并全局 atom map。
  - `build_hop_map`、`build_hop_databases`：构建回退数据库。
- `pipeline.global_bonded`
  - `build_global_bonded`：合并全局 bonded 表。
- `pipeline.parameterize`
  - `build_atom_match_with_logs`：原子匹配封装。
  - `assign_bonded_params_with_logs`：多体匹配封装。
  - `parameterize_lammps`：参数化总入口。

### 4.3 Domain
- `domain.env_key_codec`：env key 规范化与拆分列编码。
- `domain.env_key_match`：原子环境特征提取。
- `domain.atom_match`：原子级匹配与 fallback。
- `domain.keymap_merge`：keymap 合并与统计。
- `domain.term_enumeration`：拓扑项枚举。
- `domain.bonded_observed`：样本多体观测映射。
- `domain.bonded_global_merge`：全局多体合并。
- `domain.bonded_match`：多体系数匹配。

### 4.4 IO
- `io.input`：结构/LAMMPS/数据库读取解析。
- `io.output`：CSV/LAMMPS 写出。
- `io.log`：构建、匹配、冲突日志。

## 5. 关键内部常量

- 工作流默认数据库目录：`USER_DEFAULT_DB_DIR = "database"`。
- 参数化原子 fallback：主查 hop2，回退 `(1, 0)`。
- env 特征 hop 深度：内部固定为 `2`。

## 6. 开发验证

代码变更后执行：

```bash
pytest -q
```

标准夹具回归：

```bash
pytest tests/test_standard_validation.py -q
```
