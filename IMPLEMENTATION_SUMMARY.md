# MacroMapFF — 已实现内容

> 本文档记录已完成实现的功能模块，与架构设计文档 `REBUILD_PLAN.md` 区分开来。

---

## 项目结构

```
MacroMapFF/
├── src/macromapff/
│   ├── __init__.py      ✅ 已实现
│   ├── db.py            ✅ 已实现
│   ├── mol.py           ✅ 已实现
│   ├── utils.py         ✅ 已实现
│   ├── encode.py        ✅ 已实现 (Phase 2)
│   ├── lmp.py           ✅ 已实现 (Phase 3)
│   ├── fallback.py      ✅ 已实现 (Phase 5)
│   ├── cli.py           ✅ 已实现 (Phase 4/5/6)
├── tests/
│   └── test_phase1.py  ✅ 已实现
├── tests/fixtures/standard/   （测试数据）
├── pyproject.toml       ✅ 已配置
└── REBUILD_PLAN.md      （架构设计）
```

---

## 架构概览

```
数据流：

  .mol/.pdb 文件
       │
       ▼
  MolReader (mol.py) ──computeHopKeys()──► encode.py
                                               │
       │                                      │
       │  hop0/hop1/hop2 env dict             │
       ▼                                      ▼
  LammpsData (dataclass) ◄───────────── generateLammps()
       │                                      │
       │                                      │
       ▼                                      ▼
  数据库 (db.py)                         .lmp 文件
  MacroMapDB
```

**核心设计**：`LammpsData` dataclass 作为 LAMMPS 文件和工具内部表示之间的桥梁。

---

## 已实现模块

### 1. `src/macromapff/db.py` — MacroMapDB

数据库 facade，基于 pickle 存储。

| 功能 | 方法 |
|------|------|
| 加载/保存 | `load()`、`save()` |
| 原子类型（hop2级） | `insertAtomType()`、`getAtomType()` |
| Hop1 映射 | `insertHop1Key()` |
| Hop0 映射 | `insertHop0Key()` |
| 键参数 | `insertBondParam()`、`lookupBondParam()` |
| 角参数 | `insertAngleParam()`、`lookupAngleParam()` |
| 二面角参数 | `insertDihedralParam()`、`lookupDihedralParam()` |
| Improper 参数 | `insertImproperParam()`、`lookupImproperParam()` |
| 合并样本 | `mergeSample()` |
| 导出 | `export()` |

**数据库结构**（pickle 文件内部）：

```python
{
    "atom_types": {},     # hop2Key → {element, hop0_key, lammps_type, mass, sigma, epsilon, source}
    "hop1_keymap": {},    # hop1Key → {hop0_key, lammps_types: []}
    "hop0_keymap": {},    # hop0Key → {lammps_types: []}
    "bond_params": {},    # (hop0KeyA, hop0KeyB) → {k, r0}
    "angle_params": {},    # (hop0KeyA, hop0KeyB, hop0KeyC) → {k, theta0}
    "dihedral_params": {}, # (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD) → {coeffs: []}
    "improper_params": {}, # (hop0KeyCenter, hop0KeyA, hop0KeyB, hop0KeyC) → {coeffs: []}
}
```

---

### 2. `src/macromapff/encode.py` — 环境键编码 (Phase 2)

底层编码函数，直接操作 RDKit 对象。

| 函数 | 说明 |
|------|------|
| `encodeAtomEnvHop0(mol, atom)` | 计算 hop0 环境字典 |
| `encodeAtomEnvHop1(mol, atom)` | 计算 hop1 环境字典 |
| `encodeAtomEnvHop2(mol, atom)` | 计算 hop2 环境字典 |
| `encodeEnvKey(envDict)` | 将环境字典编码为 SHA-256 |
| `computeHopKeys(envHop2)` | 从 hop2 字典计算三个层的 key |
| `_neighborSignature()` | 生成邻居签名 "z:bondType:charge" |
| `_bondKindList()` | 生成 S/D/T/A/U 键类型列表 |
| `_shellNeighbors(depth)` | 获取指定深度的邻居描述符 |

---

### 3. `src/macromapff/mol.py` — MolReader

分子结构解析，基于 RDKit。调用 `encode.py` 的函数。

| 功能 | 方法 |
|------|------|
| 读取原子 | `getAtoms()` → `list[dict]` |
| 读取键 | `getBonds()` → `list[dict]` |
| 读取坐标 | `getCoords()` → `dict[int, tuple[float,float,float]]` |
| Hop0/1/2 环境 | `computeHop0/1/2Env(atomIdx)` |
| 计算三层 Key | `computeHopKeys(molReader, atomIdx)` |

---

### 4. `src/macromapff/lmp.py` — LAMMPS 文件读写 (Phase 3)

**核心数据结构**：`LammpsData` dataclass

```python
@dataclass
class LammpsData:
    header_comment: str = ""
    atoms: int = 0
    bonds: int = 0
    # ...
    masses: list[tuple[int, float]] = field(default_factory=list)
    pair_coeffs: list[tuple[int, float, float]] = field(default_factory=list)
    bond_coeffs: list[tuple[int, float, float]] = field(default_factory=list)
    angle_coeffs: list[tuple[int, float, float]] = field(default_factory=list)
    dihedral_coeffs: list[tuple[int, float, float, float, float]] = field(default_factory=list)
    improper_coeffs: list[tuple[int, ...]] = field(default_factory=list)
    atom_records: list[tuple[int, int, int, float, float, float, float]] = field(default_factory=list)
    bond_records: list[tuple[int, int, int, int]] = field(default_factory=list)
    angle_records: list[tuple[int, int, int, int, int]] = field(default_factory=list)
    dihedral_records: list[tuple[int, int, int, int, int, int]] = field(default_factory=list)
    improper_records: list[tuple[int, int, int, int, int, int]] = field(default_factory=list)
```

**主要函数**：

| 函数 | 说明 |
|------|------|
| `parseLammps(path)` → `LammpsData` | 解析 LAMMPS 文件 |
| `generateLammps(data: LammpsData, outPath)` | 将 LammpsData 写入文件 |

---

### 5. `src/macromapff/fallback.py` — 三级回退 (Phase 5)

| 函数 | 说明 |
|------|------|
| `resolveAtomType(hop2Key, hop1Key, hop0Key, element, db)` | 返回 `(lammpsType, hop0Key)` |

**回退流程**：hop2 精确 → hop1 匹配 → hop0 匹配 → 元素级兜底

---

### 6. `src/macromapff/utils.py` — 工具函数

| 函数/常量 | 说明 |
|-----------|------|
| `USER_DEFAULT_DB_DIR` | 默认数据库目录 `./database` |
| `USER_DEFAULT_DB_PATH` | 默认数据库路径 `./database/db.pkl` |
| `ensureDir(path)` | 递归创建目录 |
| `setupLogging(level)` | 配置根日志 |
| `pickleLoad(path)` | 从 pickle 反序列化 |
| `pickleSave(obj, path)` | 序列化到 pickle |

---

### 7. `src/macromapff/cli.py` — CLI 入口

三个命令：

**build-db**：
```
macromapff build-db <samples_dir> [--db-dir <path>] [-v]
```

**parameterize**：
```
macromapff parameterize <mol_file> [--out <path>] [--db-path <path>] [-v]
```

**add-samples**（暂未实现）

---

### 8. `tests/test_phase1.py` — 验收测试

17 个测试用例，覆盖 MolReader 和 MacroMapDB 的基础功能。

---

## 数据流示例

```
parameterize 流程：

1. MolReader 读取目标分子
          ↓
2. 对每个原子：
   - computeHopKeys() 计算 hop0/1/2 key
   - resolveAtomType() 三级回退 → lammps_type
          ↓
3. 构建 LammpsData 对象（原子、键、参数等）
          ↓
4. generateLammps(lmpData, outPath) 一次性写出
```

---

## 技术要点

1. **LammpsData 作为统一格式**：
   - `parseLammps` 返回 `LammpsData`
   - `parameterize` 构建 `LammpsData`
   - `generateLammps` 接收 `LammpsData` 写出
   - 类型安全、IDE 友好、可验证

2. **模块职责清晰**：
   - `encode.py` — 底层编码
   - `mol.py` — 分子读取 + 调用 encode
   - `lmp.py` — LAMMPS 文件 I/O + LammpsData 数据结构
   - `fallback.py` — 回退逻辑
   - `cli.py` — 命令编排
