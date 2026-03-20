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
      pipeline/
        build_envkey_mapping.py
        build_final_keymap.py
        build_hop_keymap.py
        extract_multiatom_terms.py
        build_multiatom_master.py
        generate_lammps_data_from_mol2.py
  tests/
  pyproject.toml
  README.md
```

说明：

- 算法脚本已迁移到 `src/macromapff/pipeline/`。
- 历史可运行工程（分子、segments、outputs、legacy scripts）已整体迁入 `examples/ps_odms7poss_legacy/`。
- 根目录仅保留开源代码与工程化配置，和 example 数据完全隔离。

## 5. 迁移策略（本次执行）

本次迁移遵循“只改工程组织，不改算法逻辑”：

1. 将历史 Python 脚本从 `scripts/` 迁移到 `src/macromapff/pipeline/`。
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

```bash
MacroMapFF build-envkey --help
MacroMapFF build-final-keymap --help
MacroMapFF generate-lammps --help
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
