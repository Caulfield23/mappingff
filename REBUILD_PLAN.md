# MacroMapFF 重构计划书

> 基于任务逻辑的全新架构，不受旧代码实现方式影响。

---

## 1. 项目目标

MacroMapFF 是一个**分子力场参数化 CLI 工具**：

- **输入**：样本分子的 `.mol` 结构文件 + 对应的 LAMMPS data file (`.lmp`)
- **输出**：新分子的参数化 LAMMPS data file

### 三个命令

| 命令 | 功能 |
|------|------|
| `build-db` | 从样本文件夹构建参数数据库 |
| `add-samples` | 追加新样本到现有数据库 |
| `parameterize` | 对目标分子进行参数化，输出 LAMMPS 文件 |

---

## 2. 数据流全景图

```
┌─────────────────────────────────────────────────────────────────┐
│                         build-db                                 │
│  samples/seg1/.mol + .lmp                                       │
│                ┌──────────────────────┐                         │
│                │  解析 mol → 分子结构  │                         │
│                │  解析 lmp → LAMMPS   │                         │
│                │    参数映射关系       │                         │
│                └──────────┬───────────┘                         │
│                           │                                      │
│                ┌──────────▼───────────┐                          │
│                │   AtomEnvEncoder    │  ← 核心：计算环境特征键   │
│                │  (hop2/hop1/hop0)   │                          │
│                └──────────┬───────────┘                         │
│                           │                                      │
│                ┌──────────▼───────────┐                          │
│                │  存入 pickle 数据库  │                          │
│                └──────────────────────┘                          │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                       parameterize                               │
│  target/PS-oDMS7POSS.mol                                        │
│                ┌──────────────────────┐                          │
│                │  解析 mol → 分子结构  │                          │
│                └──────────┬───────────┘                          │
│                           │                                      │
│                ┌──────────▼───────────┐                          │
│                │   AtomEnvEncoder     │  ← 计算目标分子环境键    │
│                │  (hop2/hop1/hop0)   │                          │
│                └──────────┬───────────┘                          │
│                           │                                      │
│                ┌──────────▼───────────┐                          │
│                │   DB Lookup +       │  ← hop2 → hop1 → hop0    │
│                │   Hop Fallback       │     三级回退匹配          │
│                └──────────┬───────────┘                          │
│                           │                                      │
│                ┌──────────▼───────────┐                          │
│                │  生成 LAMMPS data   │                          │
│                │    file (.lmp)     │                          │
│                └──────────────────────┘                          │
└─────────────────────────────────────────────────────────────────┘
```

---

## 3. 数据库设计（Pickle 方案）

### 3.1 为什么用 Pickle 而不是 SQLite

- 数据规模：上千条 atom type 记录，内存完全放得下
- 查询模式：**随机查找**（通过 env_key 查），不是范围扫描
- 增量更新：`add-samples` 需要高效插入
- 序列化/反序列化：pickle 一行代码，比 SQLite 的 C 扩展更简单

### 3.2 数据库结构

```python
# db.pkl 结构
{
    "version": 1,
    "meta": {
        "built_at": "2026-04-09T10:00:00",
        "sample_count": 4,
        "source_segments": ["segment1", "segment2", "segment3", "segment4"],
    },

    # ── 原子类型表（hop2 级别，唯一键 = hop2_key）────────────────
    # 用于 atom typing 的精确匹配表
    # key = hop2_key (SHA-256)
    "atom_types": {
        "a3f8c2d1...": {
            "element": "C",           # 元素符号
            "hop0_key": "c1f9a...",  # hop0_key（SHA-256），用于成键参数查找
            "lammps_type": 7,         # 该 hop2 类型对应的 LAMMPS atom type ID
            "mass": 12.011,
            "sigma": 3.5,              # Pair potential 参数
            "epsilon": 0.066,
            "source": ["segment1_25", "segment2_9"],  # 来源样本
        },
        # ...更多条目
    },

    # ── Hop1 KeyMap（hop1 级别，降级匹配）─────────────────────────
    # 当 hop2_key 不在数据库时，用 hop1_key 降级匹配
    "hop1_keymap": {
        # key = hop1_key (SHA-256)
        "b7d2e...": {
            "hop0_key": "c1f9a...",       # hop0_key
            "lammps_types": [3, 5, 7],     # 该 hop1 下所有 hop2 对应的 lammps_type，取最小
        },
        # ...更多条目
    },

    # ── Hop0 KeyMap（hop0 级别，再降一级）────────────────────────
    # 当 hop1_key 不在数据库时，用 hop0_key 再降一级
    "hop0_keymap": {
        # key = hop0_key (SHA-256)
        "c1f9a...": {
            "lammps_types": [3, 5],        # 该 hop0 下所有 lammps_type，取最小
        },
        # ...更多条目
    },

    # ── 成键参数表（hop0 级别匹配）──────────────────────────────
    # key = (hop0_key_a, hop0_key_b)，按字典序排列
    "bond_params": {
        ("c1f9a...", "d2e8b..."): {"k": 340.0, "r0": 1.09},
        ("c1f9a...", "e3f9c..."): {"k": 320.0, "r0": 1.41},
        # ...更多条目
    },

    # 角度参数（angle: 3个原子位置，中间固定，两端可交换）
    "angle_params": {
        ("c1f9a...", "d2e8b...", "e3f9c..."): {"k": 37.5, "theta0": 110.7},
        # ...更多条目
    },

    # 二面角参数（dihedral: 4个原子位置，两端可交换）
    "dihedral_params": {
        ("c1f9a...", "d2e8b...", "e3f9c...", "f4a0d..."): {"coeffs": [0.0, 0.0, 0.3, 0.0]},
        # ...更多条目
    },

    #  Improper 参数（4个原子，第一个是中心原子，后三个可交换）
    "improper_params": {
        ("c1f9a...", "d2e8b...", "e3f9c...", "f4a0d..."): {"coeffs": [0.0, -1, 2]},
        # ...更多条目
    },
}
```

### 3.3 关键设计决策

**hop0_key 直接作为成键参数表的 key**：
- 每个原子只有一个 hop0_key（SHA-256），不是集合
- C-O 和 O-C 在数据库中存储为 `(hop0_key_C, hop0_key_O)`，按字典序排列
- parameterize 时按相同方式构建 key 查找即可

### 3.4 为什么 env_key 用 SHA-256

- 固定 32 字节，字典查找极快
- 规范化：排序所有键保证一致性

```python
import hashlib
import json

def encodeEnvKey(envDict: dict) -> str:
    """将环境描述字典编码为固定长度哈希字符串"""
    normalized = json.dumps(envDict, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(normalized.encode()).hexdigest()
```

---

## 4. 核心算法：环境特征键编码

### 4.1 什么是"环境"

每个原子不是孤立的，它的化学环境由以下特征描述：

```
目标原子 C (z=6, degree=4)
├── 1步邻居 (hop1_shell): 3个C, 1个H
│   ├── C(sp3) - 键类型S, 无芳香性
│   ├── C(sp3) - 键类型S, 无芳香性
│   └── H      - 键类型S
├── 2步邻居 (hop2_shell): 通过hop1延伸出去的原子
└── 局部特征: 价态、是否在环中、杂化方式、电荷等
```

### 4.2 三层环境键

| 层 | 含义 | 用途 |
|----|------|------|
| `hop2` | 完整环境（hop1 + hop2 邻居） | atom typing 精确匹配，键 = hop2_key |
| `hop1` | 简化环境（只 hop1 邻居） | Fallback 中间层，键 = hop1_key |
| `hop0` | 粗粒度环境（邻居的原子序数+键类型，不含 hop2） | **成键参数匹配**，键 = hop0_key |

### 4.3 回退匹配流程

**目标**：每个原子必须分配到一个 lammps_type（用于 mass/sigma/epsilon），同时拿到 hop0_key（用于成键参数查找）。

```
parameterize 中对每个原子：
    hop2Key = encodeEnvKey(hop2Env)
    hop1Key = encodeEnvKey(hop1Env)
    hop0Key = encodeEnvKey(hop0Env)

    if hop2Key in db.atomTypes:
        lammpsType = db.atomTypes[hop2Key]["lammps_type"]
        hop0Key = db.atomTypes[hop2Key]["hop0_key"]

    elif hop1Key in db.hop1Keymap:
        lammpsType = min(db.hop1Keymap[hop1Key]["lammps_types"])
        hop0Key = db.hop1Keymap[hop1Key]["hop0_key"]

    elif hop0Key in db.hop0Keymap:
        lammpsType = min(db.hop0Keymap[hop0Key]["lammps_types"])

    else:
        # 元素级兜底：找该元素最常见的 lammps_type
        element = atom.element
        lammpsType = _fallbackByElement(element, db)
```

### 4.4 旧方案的问题

- hop_env 分成 3 个 CSV 文件（hop0/hop1/hop2_KeyMap.csv），parameterize 时要遍历
- **新方案**：所有 fallback 存在一个 dict 里，一次查找搞定
- hop0 和 hop2 的关系没有显式建模，导致成键参数匹配逻辑分散

---

## 5. 成键参数的存储与查找

### 5.1 存储（build-db 阶段）

对每个 bond，根据原子索引找到 hop0_key，按字典序排列后存储：

```python
hop0KeyA = atomHop0Key[bond["a1"]]   # "c1f9a..."
hop0KeyB = atomHop0Key[bond["a2"]]   # "d2e8b..."

# 按字典序排列，保证 C-O 和 O-C 存为相同的 key
if hop0KeyA <= hop0KeyB:
    key = (hop0KeyA, hop0KeyB)
else:
    key = (hop0KeyB, hop0KeyA)

params = bondCoeffs[bond["type_id"]]

if key in bondParams:
    # 均值合并
    bondParams[key]["k"] = (bondParams[key]["k"] + params["k"]) / 2
    bondParams[key]["r0"] = (bondParams[key]["r0"] + params["r0"]) / 2
else:
    bondParams[key] = params.copy()
```

同样逻辑适用于 angle / dihedral / improper 的存储。

### 5.2 查找（parameterize 阶段）

直接用 hop0_key 查找：

```python
def lookupBondParam(hop0KeyA: str, hop0KeyB: str, db) -> dict | None:
    if hop0KeyA <= hop0KeyB:
        key = (hop0KeyA, hop0KeyB)
    else:
        key = (hop0KeyB, hop0KeyA)
    return db.bondParams.get(key)
```

同样逻辑适用于 angle / dihedral / improper 的查找。

---

## 6. 文件结构（新架构）

```
MacroMapFF/
├── src/macromapff/
│   ├── __init__.py
│   ├── cli.py                 # 入口：argparse 分发
│   ├── db.py                  # 数据库读写（pickle）
│   ├── mol.py                 # mol 文件解析（RDKit）
│   ├── lmp.py                 # LAMMPS data file 解析/生成
│   ├── encode.py              # 环境键编码（核心算法）
│   ├── fallback.py            # Hop fallback 逻辑
│   └── utils.py               # 日志、路径等辅助
├── tests/
│   └── fixtures/standard/     # 标准测试数据（不动）
│       ├── segdata/           # 4个样本分子
│       └── target/            # 目标分子
├── pyproject.toml
└── REBUILD_PLAN.md
```

**文件数量：8 个 Python 文件（vs 旧架构 16 个）**

---

## 7. 各模块职责

### 7.1 `cli.py` — 入口

```python
# 无状态，所有逻辑委托给 db.py / mol.py / lmp.py
def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    sub.add_parser("build-db", ...)
    sub.add_parser("add-samples", ...)
    sub.add_parser("parameterize", ...)
    args = parser.parse_args()
    # 分发
```

### 7.2 `db.py` — 数据库

```python
class MacroMapDB:
    """数据库 facade，提供所有数据库操作"""
    def __init__(self, path: Path): self._path = path; self._data = None
    def load(self) -> None: ...
    def save(self) -> None: ...

    # atom_types 操作（hop2 级别）
    def insertAtomType(self, hop2Key, info) -> None: ...
    def getAtomType(self, hop2Key) -> dict | None: ...

    # hop1_keymap 操作
    def insertHop1Key(self, hop1Key, hop0Key, lammpsType) -> None: ...

    # hop0_keymap 操作
    def insertHop0Key(self, hop0Key, lammpsType) -> None: ...

    # 成键参数操作（hop0 级别）
    def insertBondParam(self, key: tuple, params: dict) -> None: ...
    def lookupBondParam(self, hop0KeyA, hop0KeyB) -> dict | None: ...
    def insertAngleParam(self, key: tuple, params: dict) -> None: ...
    def lookupAngleParam(self, ...) -> dict | None: ...
    def insertDihedralParam(self, key: tuple, params: dict) -> None: ...
    def lookupDihedralParam(self, ...) -> dict | None: ...
    def insertImproperParam(self, key: tuple, params: dict) -> None: ...
    def lookupImproperParam(self, ...) -> dict | None: ...

    def mergeSample(self, sampleData: dict) -> None: ...
    def export(self) -> dict: ...
```

### 7.3 `mol.py` — 分子结构解析

```python
from rdkit import Chem

class MolReader:
    """解析 .mol/.pdb 文件，提取分子图"""
    def __init__(self, path: Path): ...
    def getAtoms(self) -> list[dict]: ...
    def getBonds(self) -> list[dict]: ...
    def computeHop0Env(self, atomIdx) -> dict: ...  # 计算 hop0 级别环境
    def computeHop1Env(self, atomIdx) -> dict: ...  # 计算 hop1 级别环境
    def computeHop2Env(self, atomIdx) -> dict: ...   # 计算 hop2 级别环境
```

### 7.4 `encode.py` — 环境键编码（核心）

```python
def encodeAtomEnvHop0(rdkitAtom, mol) -> dict:
    """hop0 级别环境：中心原子 + hop1 邻居（不含 hop2）"""
    return {
        "z": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "degree": atom.GetTotalDegree(),
        "hybridization": str(atom.GetHybridization()),
        "in_ring": int(atom.IsInRing()),
        "ring_count": atom.GetRingInfo().NumRings(),
        "total_hs": atom.GetTotalNumHs(),
        "neighbor_sig": _neighborSignature(mol, atom),
        "bond_kinds": _bondKindList(mol, atom),
    }

def encodeAtomEnvHop1(rdkitAtom, mol) -> dict:
    """hop1 级别环境：hop0 + hop1_shell 详细信息"""
    env = encodeAtomEnvHop0(rdkitAtom, mol)
    env["hop1_shell"] = _shellNeighbors(mol, atom, depth=1)
    return env

def encodeAtomEnvHop2(rdkitAtom, mol) -> dict:
    """hop2 级别环境：hop1 + hop2_shell 详细信息"""
    env = encodeAtomEnvHop1(rdkitAtom, mol)
    env["hop2_shell"] = _shellNeighbors(mol, atom, depth=2)
    return env

def computeHopKeys(envHop2: dict) -> tuple[str, str, str]:
    """从 hop2 环境生成 hop2/hop1/hop0 三层 SHA-256 key"""
    hop2Key = sha256(json.dumps(envHop2, sort_keys=True))
    hop1Env = _stripHop2(envHop2)
    hop1Key = sha256(json.dumps(hop1Env, sort_keys=True))
    hop0Env = _stripHop1(envHop1)
    hop0Key = sha256(json.dumps(hop0Env, sort_keys=True))
    return hop2Key, hop1Key, hop0Key
```

### 7.5 `lmp.py` — LAMMPS 文件

```python
def parseLammps(path: Path) -> dict:
    """解析 LAMMPS data file"""
    return {
        "atoms": atomRecords,              # [{idx, type_id, charge, x, y, z}]
        "bonds": bondRecords,
        "angles": angleRecords,
        "dihedrals": dihedralRecords,
        "impropers": improperRecords,
        "pairCoeffs": pairCoeffs,
        "bondCoeffs": bondCoeffs,
        "angleCoeffs": angleCoeffs,
        "dihedralCoeffs": dihedralCoeffs,
        "improperCoeffs": improperCoeffs,
    }

def generateLammps(molPath: Path, db: MacroMapDB, outPath: Path,
                    chargeStrategy: ChargeStrategy | None = None) -> None:
    """生成参数化的 LAMMPS data file

    Args:
        chargeStrategy: 电荷处理策略接口，默认不使用样本电荷（留空，后续实现）
    """
    # 1. 读 mol
    # 2. 对每个原子编码 + 查 db
    # 3. 分配 lammps_type（连续编号）
    # 4. 用 hop0_key 查成键参数（bond/angle/dihedral/improper）
    # 5. 组装输出（电荷字段由 chargeStrategy 填充，暂不处理）
```

**电荷接口**（后续实现）：
```python
class ChargeStrategy(Protocol):
    """电荷策略接口，由用户后续实现"""
    def computeCharges(self, molPath: Path, atomTypeMap: dict) -> dict[int, float]:
        """返回 atom_idx → charge 的映射"""
        ...
```

### 7.6 `utils.py` — 工具函数

日志配置、默认路径（USER_DEFAULT_DB_DIR）、pickle 序列化辅助函数等。

```python
USER_DEFAULT_DB_DIR = Path("./database")

def ensureDir(path: Path) -> None: ...

def setupLogging(level: int = logging.INFO) -> None: ...
```

### 7.7 `fallback.py` — 回退逻辑

```python
def resolveAtomType(hop2Key: str, hop1Key: str, hop0Key: str,
                      element: str, db: MacroMapDB) -> tuple[int, str]:
    """解析原子类型，返回 (lammpsType, hop0Key)

    依次尝试 hop2 → hop1 → hop0 → 元素级兜底
    """
    if hop2Key in db.atomTypes:
        at = db.atomTypes[hop2Key]
        return at["lammps_type"], at["hop0_key"]

    if hop1Key in db.hop1Keymap:
        entry = db.hop1Keymap[hop1Key]
        return min(entry["lammps_types"]), entry["hop0_key"]

    if hop0Key in db.hop0Keymap:
        return min(db.hop0Keymap[hop0Key]["lammps_types"]), hop0Key

    # 最终兜底：元素级
    return _fallbackByElement(element, db)
```

---

## 8. 工作流详解

### 8.1 `build-db` 流程

```
samples/
├── segment1/
│   ├── segment1.mol
│   └── segment1.lammps.lmp
└── ...

Step 1: 遍历每个 segment 文件夹
        molReader = MolReader(segment.mol)
        lmpData = parseLammps(segment.lammps.lmp)

Step 2: 对每个原子（共 N 个）
        # 计算三层环境键
        envHop0 = molReader.computeHop0Env(atomIdx)
        envHop1 = molReader.computeHop1Env(atomIdx)  # = envHop0 + hop1_shell
        envHop2 = molReader.computeHop2Env(atomIdx)  # = envHop1 + hop2_shell
        hop2Key, hop1Key, hop0Key = computeHopKeys(envHop2)
        atomType = lmpData["atom_type_map"][atomIdx]

        # 存储到 atom_types（hop2 级别）
        db.insertAtomType(hop2Key, {
            "element": element,
            "hop0_key": hop0Key,
            "lammps_type": atomType,
            "mass": lmpData["pair_coeffs"][atomType]["mass"],
            "sigma": lmpData["pair_coeffs"][atomType]["sigma"],
            "epsilon": lmpData["pair_coeffs"][atomType]["epsilon"],
            "source": [f"{segName}_{atomIdx}"]
        })

        # 存储到 hop1_keymap
        db.insertHop1Key(hop1Key, hop0Key, atomType)

        # 存储到 hop0_keymap
        db.insertHop0Key(hop0Key, atomType)

Step 3: 对每个成键参数（bond / angle / dihedral / improper）
        从 lmpData 提取原子上的键参数

        例：bond C(hop0Key="c1f9a...") - H(hop0Key="d2e8b...")
        → 按字典序排列 → ("c1f9a...", "d2e8b...")
        → 存储 k, r0

Step 4: 合并同键（同一 hop2Key 出现多次时更新 source 和均值）
        if hop2Key in db.atomTypes:
            existing = db.atomTypes[hop2Key]
            # 更新 sigma/epsilon/mass 均值
            # 追加到 source

        # hop1_keymap 和 hop0_keymap 中同一 key 有多个 lammpsType 时取最小
        # 成键参数同 key 时取均值

Step 5: pickle.dump(db.export(), open(dbPath, "wb"))
```

### 8.2 `add-samples` 流程

```
Step 1: 加载现有 db.pkl
Step 2: 解析新样本（同 build-db Step 2-3）
Step 3: 增量合并
        for hop2Key, info in newAtomTypes:
            if hop2Key in existing:
                _mergeAtomType(existing[hop2Key], info)
            else:
                existing[hop2Key] = info

        # 成键参数：同 key → 取参数均值

Step 4: pickle.dump(...)
```

### 8.3 `parameterize` 流程

```
target/PS-oDMS7POSS.mol

Step 1: 加载 db.pkl（整个数据库进内存）
Step 2: 解析目标分子 mol 文件
Step 3: 对每个原子
        envHop2 = molReader.computeHop2Env(atomIdx)
        hop2Key, hop1Key, hop0Key = computeHopKeys(envHop2)
        lammpsType, hop0Key = resolveAtomType(hop2Key, hop1Key, hop0Key, element, db)
        → 记录 atomIdx → (lammpsType, hop0Key)

Step 4: 分配新的 LAMMPS type ID
        uniqueTypes = sorted(set(lammpsTypes))
        typeMapping = {old: new for new, old in enumerate(uniqueTypes, 1)}

Step 5: 查找成键参数（直接用 hop0Key 匹配）
        for bond in targetBonds:
            hop0KeyA = atomHop0Key[bond.a1]
            hop0KeyB = atomHop0Key[bond.a2]
            params = db.lookupBondParam(hop0KeyA, hop0KeyB)
            # 若找不到 → 警告

        # angle, dihedral, improper 同理

Step 6: 生成 LAMMPS data file
        writeHeader(atomCount, bondCount, ...)
        writeMasses(uniqueTypes, db)
        writePairCoeffs(uniqueTypes, db)
        writeBonds(bondCount, bondParams, ...)
        writeAngles(...)
        writeDihedrals(...)
        writeImpropers(...)
        # 电荷字段由 chargeStrategy 填充（暂不处理）
```

---

## 9. 实现顺序

### Phase 1：基础设施

```
1. pyproject.toml 配置
   - rdkit-python 依赖

2. db.py — 数据库读写
   - 纯 pickle，无复杂逻辑
   - load() / save() / insert_* / get_*

3. mol.py — 分子读取
   - 用 RDKit 解析 mol/pdb 文件
   - 提取原子、键、坐标、compute_hop0/1/2_env
```

**验收标准**：能够读取一个 .mol 文件，拿到原子列表、键列表，以及任意原子的 hop0/hop1/hop2 环境。

### Phase 2：核心算法

```
4. encode.py — 环境键编码
   - 实现 compute_hop_keys()
   - 目标：对 segment2.mol 中的原子，能生成 hop2/hop1/hop0 三层 key

```

### Phase 3：LAMMPS 文件

```
6. lmp.py — LAMMPS 文件读写
   - 解析 .lmp 文件，提取 atom type / mass / sigma/epsilon / bond/angle/dihedral/improper 参数
   - 生成新的 LAMMPS data file
```

**验收标准**：
- `parseLammps(segment2.lammps.lmp)` 能正确提取所有参数
- `generateLammps()` 输出的文件格式与参考一致

### Phase 4：build-db

```
7. build-db 命令实现
   - 遍历 samples/ 下所有 segment
   - 对每个 segment 解析 mol + lmp
   - 用 encode.py 计算环境键
   - 用 db.py 存储 atomTypes、hop1Keymap、hop0Keymap
   - 成键参数按 hop0Key 字典序排列后存储
```

**验收标准**：`macromapff build-db tests/fixtures/standard/segdata` 能生成 db.pkl，且：
- atom_types 包含所有 hop2 键
- hop1_keymap 和 hop0_keymap 包含降级路径
- bond_params 数量合理（合并后比原始样本少）

### Phase 5：parameterize + fallback

```
8. fallback.py — 回退逻辑
   - 实现 resolveAtomType
   - hop2 → hop1 → hop0 三级回退

9. parameterize 命令实现
   - 加载 db.pkl
   - 解析目标分子
   - 对每个原子做环境编码 + 三级回退查找
   - 用 hop0Key 查找成键参数
   - 生成 LAMMPS data file
```

**验收标准**：`macromapff parameterize tests/fixtures/standard/target/PS-oDMS7POSS.mol` 生成 `PS-oDMS7POSS_param.lmp`，与 tests/artifacts/standard/case_standard_workflow/PS-oDMS7POSS_param.lmp 逐原子类型对比一致。

### Phase 6：add-samples

```
10. add-samples 命令实现
    - 加载现有 db.pkl
    - 增量插入新样本
    - pickle.dump
```

### Phase 7：测试

```
11. 写测试用例
    - testBuildDb: 验证数据库正确性
    - testParameterize: 对比输出与标准 artifacts
    - testAddSamples: 验证增量插入正确
```

---

## 10. 测试策略

### 标准测试数据（不动）

```
tests/fixtures/standard/
├── segdata/
│   ├── segment1/ (.mol + .lmp)
│   ├── segment2/ (.mol + .lmp)
│   ├── segment3/ (.mol + .lmp)
│   └── segment4/ (.mol + .lmp)
└── target/
    └── PS-oDMS7POSS.mol
```

**参考输出**（对比基准）：
```
tests/artifacts/standard/case_standard_workflow/
├── PS-oDMS7POSS_param.lmp
└── atom_index_key_types.csv
```

### 测试用例

```python
def testBuildDbStandard(tmp_path):
    """build-db 能正确处理 4 个 segment 并合并"""
    dbPath = tmpPath / "db.pkl"
    result = buildDb(Path("tests/fixtures/standard/segdata"), dbPath)
    assert result["samplesCount"] == 4

    db = MacroMapDB(dbPath)
    db.load()
    assert len(db.atomTypes) > 0
    assert len(db.hop1Keymap) > 0
    assert len(db.hop0Keymap) > 0
    assert len(db.bondParams) > 0


def testParameterizeStandard(tmp_path):
    """parameterize 输出与已知正确结果一致"""
    dbPath = tmpPath / "db.pkl"
    buildDb(Path("tests/fixtures/standard/segdata"), dbPath)

    outPath = tmpPath / "out.lmp"
    parameterize(Path("tests/fixtures/standard/target/PS-oDMS7POSS.mol"), dbPath, outPath)

    expected = Path("tests/artifacts/standard/case_standard_workflow/PS-oDMS7POSS_param.lmp")
    assert outPath.readText() == expected.readText()


def testAddSamplesIncremental(tmp_path):
    """add-samples 只追加新类型"""
    dbPath = tmpPath / "db.pkl"
    buildDb(Path("segdata/segment1-3"), dbPath)
    countBefore = len(MacroMapDB(dbPath).load().atomTypes)

    addSamples(Path("segdata/segment4"), dbPath)
    db = MacroMapDB(dbPath).load()
    assert len(db.atomTypes) >= countBefore
```

---

## 11. 关键设计决策

### 11.1 为什么 env_key 用 SHA-256

- **可哈希**：dict 不能作为 dict 的 key，JSON 字符串可以
- **长度固定**：32 字节 vs 几十到几百字节的 JSON 字符串
- **二进制比较更快**：不需要解析 JSON

### 11.2 为什么不用 SQLite

| | Pickle + dict | SQLite |
|--|--|--|
| 读取 | 全量加载 O(n) | 按查询 O(log n) |
| 写入 | 整体重写 O(n) | 增量插入 O(1) |
| 依赖 | Python 标准库 | 需要 C 扩展 |
| 适合规模 | ~10k 条记录 | ~100k+ 条记录 |

对于这个工具的规模（上千条 atom types），pickle + 全量内存加载完全够用，且更简单。

### 11.3 为什么成键参数在 hop0 级别而不是 hop2 级别

因为同一个 C-H 键在数据库中可能对应多个 hop2 类型（不同 C 的 hop2 环境），但它们都应该匹配到同一个 C-H bond 参数。在 hop2 级别存储会导致参数表膨胀且存在大量冗余。

### 11.4 电荷处理

**当前设计**：样本中的电荷**不写入数据库**，也不用于 parameterize 输出。

电荷是分子整体的属性，不能直接从样本迁移到目标分子。后续由用户通过 `ChargeStrategy` 接口实现自定义电荷计算逻辑。

```python
# 生成 LAMMPS 时，电荷部分由 charge_strategy 填充（暂不实现）
# pair_coeffs 中的 sigma/epsilon 来自数据库
# mass 来自数据库
```

### 11.5 成键参数合并策略

当同一成键模式在多个样本中出现时，取参数**均值**：
- Bond: `k` 和 `r0` 分别取均值
- Angle: `k` 和 `theta0` 分别取均值
- Dihedral/Improper: 各系数分别取均值

---

## 12. 项目结构（最终）

```
src/macromapff/
├── __init__.py
├── cli.py          # 入口
├── db.py           # pickle 数据库
├── mol.py          # RDKit mol 解析
├── lmp.py          # LAMMPS 文件读写
├── encode.py       # 环境键编码
├── fallback.py      # 三级回退解析
└── utils.py        # 日志、路径等辅助
```

### 依赖

```toml
[project]
requires-python = ">=3.10"
dependencies = [
    "rdkit>=2024.0.0",
]
```

### CLI 命令

```bash
macromapff build-db <samplesDir> [--db-dir <path>]
macromapff add-samples <samplesDir> [--db-dir <path>]
macromapff parameterize <molFile> [--out <path>] [--db-dir <path>]
```