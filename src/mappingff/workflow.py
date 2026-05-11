"""Core workflow functions for mappingff.

This module provides the main workflow functions:
    - buildDb: Build parameter database from sample molecules
    - parameterize: Parameterize target molecule and generate LAMMPS file

These functions can be imported and used programmatically, or via the CLI.
"""

from __future__ import annotations

import logging
from pathlib import Path

from mappingff.db import MacroMapDB
from mappingff.encode import getHop0Subgraph, getHop3Subgraph
from mappingff.fallback import resolveAtomType
from mappingff.lmp import LammpsData, adjustTotalCharge, generateLammps, parseLammps
from mappingff.mol import MolReader, computeHopKeys


def buildDb(samplesDir: Path, dbPath: Path) -> None:
    """Build a parameter database from sample molecules.

    Iterates over subdirectories in samplesDir, each expected to contain
    a .mol (or .pdb) and a .lmp file with the same base name.

    Args:
        samplesDir: Directory containing sample subdirectories.
        dbPath: Path to output the SQLite database file.

    """
    log = logging.getLogger("build-db")

    # Find all sample directories
    sampleDirs = sorted([d for d in samplesDir.iterdir() if d.is_dir()])

    # Initialize database
    db = MacroMapDB(dbPath)
    db.load()

    totalAtoms = 0
    totalBonds = 0
    totalAngles = 0
    totalDihedrals = 0
    totalImpropers = 0

    # Collect all items to process: (lmpPath, topoPath or None)
    items: list[tuple[Path, Path | None]] = []
    for sampleDir in sampleDirs:
        topoPath = next(sampleDir.glob("*.mol"), None) or next(
            sampleDir.glob("*.pdb"), None
        )
        lmpPath = next(sampleDir.glob("*.lmp"), None)
        if topoPath and lmpPath:
            items.append((lmpPath, topoPath))
    for lmpPath in samplesDir.glob("*.lmp"):
        items.append((lmpPath, None))

    for lmpPath, topoPath in items:
        segName = lmpPath.stem

        file_part = f" + {topoPath.name}" if topoPath else " (lmp only)"
        log.info(f"Processing {segName}: {lmpPath.name}{file_part}")

        lmpData = parseLammps(lmpPath)
        molReader = MolReader(topoPath if topoPath else lmpPath)

        atoms = molReader.getAtoms()
        bonds = molReader.getBonds()
        angles = molReader.getAngles()
        dihedrals = molReader.getDihedrals()
        impropers = molReader.getImpropers()
        totalAtoms += len(atoms)
        totalBonds += len(bonds)
        totalAngles += len(angles)
        totalDihedrals += len(dihedrals)
        totalImpropers += len(impropers)

        # Build atom_idx -> hop0_key mapping for bonded parameter lookup
        atomHop0Key: dict[int, str] = {}

        # Process each atom
        for atom in atoms:
            atomIdx = atom["idx"]
            hop3Key, hop2Key, hop1Key, hop0Key = computeHopKeys(molReader, atomIdx)
            rd_atom = molReader.mol.GetAtomWithIdx(atomIdx - 1)
            hop0_graph = getHop0Subgraph(molReader.mol, rd_atom)
            hop3_graph = getHop3Subgraph(molReader.mol, rd_atom)
            atomHop0Key[atomIdx] = hop0Key

            # Get LAMMPS type from the sample data
            atomRecord = next(
                (rec for rec in lmpData.atom_records if rec[0] == atomIdx),
                None,
            )
            if atomRecord is None:
                log.warning(f"Atom {atomIdx} not found in LAMMPS data for {segName}")
                continue

            # atomRecord format: [atom_id, mol_tag, type_id, charge, x, y, z]
            lammpsType = atomRecord[2]

            # Get pair coeffs for this type
            pairCoeff = next(
                (rec for rec in lmpData.pair_coeffs if rec[0] == lammpsType),
                None,
            )
            mass = next(
                (rec for rec in lmpData.masses if rec[0] == lammpsType),
                None,
            )

            # Insert atom type at hop3 level (finest classification)
            # lammps_type is auto-generated to ensure uniqueness per hop3Key
            db.insertAtomType(
                {
                    "hop3_key": hop3Key,
                    "hop2_key": hop2Key,
                    "hop1_key": hop1Key,
                    "hop0_key": hop0Key,
                    "element": atom["symbol"],
                    "hop0_graph": hop0_graph,
                    "hop3_graph": hop3_graph,
                    "mass": mass[1] if mass else 0.0,
                    "sigma": pairCoeff[2] if pairCoeff else 0.0,
                    "epsilon": pairCoeff[1] if pairCoeff else 0.0,
                    "charge": atomRecord[3],
                    "source": f"{segName}_{atomIdx}",
                },
            )

            # Get the auto-generated lammps_type from database
            inserted_info = db.getAtomType(hop3Key)
            assigned_lammps_type = inserted_info["lammps_type"]

            # Insert hop keymap (for external validation only)
            db.insertHopKey([hop0Key, hop1Key, hop2Key], assigned_lammps_type)

        # Process bonds
        # bond_records: [bond_id, bond_type, a1, a2]
        for bondRec in lmpData.bond_records:
            bondType = bondRec[1]
            a1 = bondRec[2]
            a2 = bondRec[3]

            hop0KeyA = atomHop0Key[a1]
            hop0KeyB = atomHop0Key[a2]

            # bondCoeffs: [type_id, K, r0] (bond_style harmonic)
            bondCoeff = next(rec for rec in lmpData.bond_coeffs if rec[0] == bondType)

            # Canonical key (order-independent)
            key = min((hop0KeyA, hop0KeyB), (hop0KeyB, hop0KeyA))

            db.insertBondParam(key, {"k": bondCoeff[1], "r0": bondCoeff[2]})

        # Process angles
        # angle_records: [angle_id, angle_type, a1, a2, a3]
        for angleRec in lmpData.angle_records:
            angType = angleRec[1]
            a1 = angleRec[2]
            a2 = angleRec[3]
            a3 = angleRec[4]

            hop0KeyA = atomHop0Key[a1]
            hop0KeyB = atomHop0Key[a2]
            hop0KeyC = atomHop0Key[a3]

            # angleCoeffs: [type_id, K, theta0] (angle_style harmonic)
            angleCoeff = next(rec for rec in lmpData.angle_coeffs if rec[0] == angType)

            # Canonical key (outer atoms swappable)
            key = min((hop0KeyA, hop0KeyB, hop0KeyC), (hop0KeyC, hop0KeyB, hop0KeyA))

            db.insertAngleParam(key, {"k": angleCoeff[1], "theta0": angleCoeff[2]})

        # Process dihedrals
        # dihedral_records: [dih_id, dih_type, a1, a2, a3, a4]
        for dihRec in lmpData.dihedral_records:
            dihType = dihRec[1]
            a1 = dihRec[2]
            a2 = dihRec[3]
            a3 = dihRec[4]
            a4 = dihRec[5]

            hop0KeyA = atomHop0Key[a1]
            hop0KeyB = atomHop0Key[a2]
            hop0KeyC = atomHop0Key[a3]
            hop0KeyD = atomHop0Key[a4]

            # dihCoeffs: [type_id, K1, K2, K3, K4] (dihedral_style opls)
            dihCoeff = next(rec for rec in lmpData.dihedral_coeffs if rec[0] == dihType)

            # Canonical key (A,B,C,D) and (D,C,B,A) are equivalent, pick lexicographically smaller
            key = min(
                (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD),
                (hop0KeyD, hop0KeyC, hop0KeyB, hop0KeyA),
            )

            db.insertDihedralParam(key, {"coeffs": dihCoeff[1:]})

        # Process impropers
        # improper_records: [imp_id, imp_type, a1, a2, a3, a4]
        for impRec in lmpData.improper_records:
            impType = impRec[1]
            a1 = impRec[2]
            a2 = impRec[3]
            a3 = impRec[4]
            a4 = impRec[5]

            hop0KeyA = atomHop0Key[a1]
            hop0KeyB = atomHop0Key[a2]
            hop0KeyC = atomHop0Key[a3]
            hop0KeyD = atomHop0Key[a4]

            # impCoeffs: [type_id, K, d, n] (improper_style cvff)
            impCoeff = next(rec for rec in lmpData.improper_coeffs if rec[0] == impType)

            # Canonical key (center fixed, others sorted)
            others = sorted([hop0KeyB, hop0KeyC, hop0KeyD])
            key = (hop0KeyA, others[0], others[1], others[2])

            db.insertImproperParam(key, {"coeffs": impCoeff[1:]})

        log.info(
            f"  {segName}: {len(atoms)} atoms, {len(bonds)} bonds, "
            f"{len(angles)} angles, {len(dihedrals)} dihedrals, {len(impropers)} impropers"
        )

        # Update metadata for this sample
        db.setMeta(
            segName,
            f"{lmpPath} | {topoPath if topoPath else '(lmp only)'}",
        )

    # Update metadata
    db.setMeta("sample_count", str(len(items)))

    # Save database
    db.save()

    log.info(f"Database saved to {dbPath}")
    log.info(f"Samples: {len(items)}")
    log.info(
        f"Processed: {totalAtoms} atoms, {totalBonds} bonds, "
        f"{totalAngles} angles, {totalDihedrals} dihedrals, {totalImpropers} impropers"
    )
    log.info(f"Hop3 atom types: {len(db.atomTypes)}")
    log.info(f"Hop2 atom types: {len(db.hop2Keymap)}")
    log.info(f"Hop1 atom types: {len(db.hop1Keymap)}")
    log.info(f"Hop0 atom types: {len(db.hop0Keymap)}")
    log.info(f"Bond types: {len(db.bondParams)}")
    log.info(f"Angle types: {len(db.angleParams)}")
    log.info(f"Dihedral types: {len(db.dihedralParams)}")
    log.info(f"Improper types: {len(db.improperParams)}")


def parameterize(
    topoPath: Path,
    dbPath: Path,
    outPath: Path | None = None,
    total_charge: float = 0.0,
) -> dict:
    """Parameterize a target molecule using the database.

    Args:
        topoPath: Path to target molecule .mol or .pdb file.
        dbPath: Path to the SQLite database file.
        outPath: Output LAMMPS file path. If None, uses <mol_file>.lmp.
        total_charge: Target total charge for the system (default: 0).
            Charge will be adjusted evenly across all non-hydrogen atoms.

    Returns:
        Dictionary with parameterization statistics.
    """
    log = logging.getLogger("parameterize")

    # Load database
    db = MacroMapDB(dbPath)
    db.load()
    log.info(f"Database loaded from {dbPath}")
    log.debug(f"  Hop3 atom types: {len(db.atomTypes)}")
    log.debug(f"  Hop2 atom types: {len(db.hop2Keymap)}")
    log.debug(f"  Hop1 atom types: {len(db.hop1Keymap)}")
    log.debug(f"  Hop0 atom types: {len(db.hop0Keymap)}")
    log.debug(f"  Bond types: {len(db.bondParams)}")
    log.debug(f"  Angle types: {len(db.angleParams)}")
    log.debug(f"  Dihedral types: {len(db.dihedralParams)}")
    log.debug(f"  Improper types: {len(db.improperParams)}")

    # Parse target molecule
    molReader = MolReader(topoPath)
    atoms = molReader.getAtoms()
    bonds = molReader.getBonds()
    angles = molReader.getAngles()
    dihedrals = molReader.getDihedrals()
    impropers = molReader.getImpropers()
    coords = molReader.getCoords()

    log.info(f"Parameterizing molecule: {topoPath.name}")

    # Determine output path
    if outPath is None:
        outPath = topoPath.with_suffix(".lmp")

    # Resolve atom types
    atomTypeMap: dict[int, int] = {}  # atomIdx -> lammps_type (resolved from db)
    atomTypeParams: dict[int, dict] = (
        {}
    )  # lammps_type -> {type: consecutive_output_type, element, mass, sigma, epsilon, charge}
    atomHop0Key: dict[int, str] = {}  # atomIdx -> hop0Key (resolved from db)
    atomHop3Key: dict[int, str] = {}  # atomIdx -> hop3Key (for debugging)
    atomFallbackLevel: dict[int, str] = {}  # atomIdx -> fallback level (hop2/hop1/hop0)
    nextAtomTypeId = 1
    hop3Matches = 0
    hop2Matches = 0
    hop1Matches = 0
    hop0Matches = 0
    noMatch = 0

    for atom in atoms:
        atomIdx = atom["idx"]
        hop3Key, hop2Key, hop1Key, hop0Key = computeHopKeys(molReader, atomIdx)
        atomHop0Key[atomIdx] = hop0Key
        atomHop3Key[atomIdx] = hop3Key

        lammpsType, resolvedHop0Key, hopLevel, matchedHop3Key = resolveAtomType(
            hop3Key, hop2Key, hop1Key, hop0Key, db
        )
        # Use lammpsType as the type identifier
        atomTypeMap[atomIdx] = lammpsType
        atomHop0Key[atomIdx] = resolvedHop0Key

        # Collect atom type parameters (deduplicated by db_lammps_type)
        if lammpsType is not None and lammpsType not in atomTypeParams:
            info = db.atomTypes[matchedHop3Key]
            atomTypeParams[lammpsType] = {
                "type": nextAtomTypeId,
                "element": info["element"],
                "mass": info["mass"],
                "sigma": info["sigma"],
                "epsilon": info["epsilon"],
                "charge": info["charge"],
            }
            nextAtomTypeId += 1

        # Track match statistics
        if hopLevel == "hop3":
            hop3Matches += 1
            log.debug(
                f"  Atom {atomIdx}: element={atom['symbol']}, lammps_type={lammpsType}"
            )
        elif hopLevel == "hop2":
            hop2Matches += 1
            atomFallbackLevel[atomIdx] = "hop2"
            log.debug(
                f"  Atom {atomIdx}: element={atom['symbol']}, lammps_type={lammpsType} (hop2 fallback)"
            )
        elif hopLevel == "hop1":
            hop1Matches += 1
            atomFallbackLevel[atomIdx] = "hop1"
            log.debug(
                f"  Atom {atomIdx}: element={atom['symbol']}, lammps_type={lammpsType} (hop1 fallback)"
            )
        elif hopLevel == "hop0":
            hop0Matches += 1
            atomFallbackLevel[atomIdx] = "hop0"
            log.debug(
                f"  Atom {atomIdx}: element={atom['symbol']}, lammps_type={lammpsType} (hop0 fallback)"
            )
        else:
            noMatch += 1
            atomFallbackLevel[atomIdx] = "none"
            log.debug(
                f"  Atom {atomIdx}: element={atom['symbol']}, lammps_type={lammpsType} (NO MATCH)"
            )

    # Look up bond parameters
    bondTypeMap: dict[int, int] = {}  # bondIdx -> newBondTypeId
    bondTypeParams: dict[int, tuple] = {}  # newBondTypeId -> (k, r0)
    nextBondTypeId = 1
    bondNoMatch = 0
    bondNoMatchAtoms: list[tuple[int, int]] = []

    for bond in bonds:
        bondIdx = bond["idx"]
        a1 = bond["a1"]
        a2 = bond["a2"]
        hop0KeyA = atomHop0Key.get(a1)
        hop0KeyB = atomHop0Key.get(a2)
        if hop0KeyA is None or hop0KeyB is None:
            raise ValueError(
                f"Bond {bondIdx} references missing atom {a1 if hop0KeyA is None else a2}"
            )

        # Look up bond param
        bondParam = db.lookupBondParam(hop0KeyA, hop0KeyB)
        if bondParam is None:
            bondNoMatch += 1
            bondNoMatchAtoms.append((a1, a2))
            continue

        # Create bond type if not seen
        paramKey = (bondParam["k"], bondParam["r0"])
        bondTypeId = None
        for btid, (bk, br0) in bondTypeParams.items():
            if abs(bk - paramKey[0]) < 0.01 and abs(br0 - paramKey[1]) < 0.001:
                bondTypeId = btid
                break

        if bondTypeId is None:
            bondTypeId = nextBondTypeId
            nextBondTypeId += 1
            bondTypeParams[bondTypeId] = paramKey

        bondTypeMap[bondIdx] = bondTypeId

    # Look up angle parameters
    angleTypeMap: dict[int, int] = {}  # angleIdx -> newAngleTypeId
    angleTypeParams: dict[int, tuple] = {}  # newAngleTypeId -> (k, theta0)
    nextAngleTypeId = 1
    angleNoMatch = 0
    angleNoMatchAtoms: list[tuple[int, int, int]] = []

    for angle in angles:
        angleIdx = angle["idx"]
        a1 = angle["a1"]
        a2 = angle["a2"]
        a3 = angle["a3"]
        hop0KeyA = atomHop0Key.get(a1)
        hop0KeyB = atomHop0Key.get(a2)
        hop0KeyC = atomHop0Key.get(a3)
        if hop0KeyA is None or hop0KeyB is None or hop0KeyC is None:
            raise ValueError(
                f"Angle {angleIdx} references missing atom {a1 if hop0KeyA is None else a2 if hop0KeyB is None else a3}"
            )

        # Look up angle param
        angleParam = db.lookupAngleParam(hop0KeyA, hop0KeyB, hop0KeyC)
        if angleParam is None:
            angleNoMatch += 1
            angleNoMatchAtoms.append((a1, a2, a3))
            continue

        # Create angle type if not seen
        paramKey = (angleParam["k"], angleParam["theta0"])
        angleTypeId = None
        for atid, (ak, ath) in angleTypeParams.items():
            if abs(ak - paramKey[0]) < 0.01 and abs(ath - paramKey[1]) < 0.1:
                angleTypeId = atid
                break

        if angleTypeId is None:
            angleTypeId = nextAngleTypeId
            nextAngleTypeId += 1
            angleTypeParams[angleTypeId] = paramKey

        angleTypeMap[angleIdx] = angleTypeId

    # Look up dihedral parameters
    dihedralTypeMap: dict[int, int] = {}  # dihIdx -> newDihTypeId
    dihedralTypeParams: dict[int, tuple] = {}  # newDihTypeId -> coeffs tuple
    nextDihedralTypeId = 1
    dihedralNoMatch = 0
    dihedralNoMatchAtoms: list[tuple[int, int, int, int]] = []

    for dihedral in dihedrals:
        dihIdx = dihedral["idx"]
        a1 = dihedral["a1"]
        a2 = dihedral["a2"]
        a3 = dihedral["a3"]
        a4 = dihedral["a4"]
        hop0KeyA = atomHop0Key.get(a1)
        hop0KeyB = atomHop0Key.get(a2)
        hop0KeyC = atomHop0Key.get(a3)
        hop0KeyD = atomHop0Key.get(a4)
        if hop0KeyA is None or hop0KeyB is None or hop0KeyC is None or hop0KeyD is None:
            raise ValueError(
                f"Dihedral {dihIdx} references missing atom {a1 if hop0KeyA is None else a2 if hop0KeyB is None else a3 if hop0KeyC is None else a4}"
            )

        # Look up dihedral param
        dihedralParam = db.lookupDihedralParam(hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD)
        if dihedralParam is None:
            dihedralNoMatch += 1
            dihedralNoMatchAtoms.append((a1, a2, a3, a4))
            continue

        # Create dihedral type if not seen
        coeffs = tuple(dihedralParam.get("coeffs", []))
        dihedralTypeId = None
        for dtid, existing_coeffs in dihedralTypeParams.items():
            if len(coeffs) == len(existing_coeffs):
                match = all(abs(c - e) < 0.01 for c, e in zip(coeffs, existing_coeffs))
                if match:
                    dihedralTypeId = dtid
                    break

        if dihedralTypeId is None:
            dihedralTypeId = nextDihedralTypeId
            nextDihedralTypeId += 1
            dihedralTypeParams[dihedralTypeId] = coeffs

        dihedralTypeMap[dihIdx] = dihedralTypeId

    # Look up improper parameters
    improperTypeMap: dict[int, int] = {}  # impIdx -> newImpTypeId
    improperTypeParams: dict[int, tuple] = {}  # newImpTypeId -> coeffs tuple
    nextImproperTypeId = 1
    improperNoMatch = 0
    improperNoMatchAtoms: list[tuple[int, int, int, int]] = []

    # Build atom and atom degree index for cvff improper filtering
    atom_name_by_idx = {a["idx"]: a for a in atoms}
    atom_degree_by_idx = {a["idx"]: a["degree"] for a in atoms}

    for improper in impropers:
        impIdx = improper["idx"]
        a1 = improper["a1"]
        a2 = improper["a2"]
        a3 = improper["a3"]
        a4 = improper["a4"]

        # improper filter: center atom (a1) must be C or N with exactly 3 neighbors
        center_atom = atom_name_by_idx.get(a1)
        if center_atom is None:
            continue
        if center_atom["symbol"] not in ("C", "N"):
            continue
        if atom_degree_by_idx.get(a1, 0) != 3:
            continue

        hop0KeyA = atomHop0Key.get(a1)
        hop0KeyB = atomHop0Key.get(a2)
        hop0KeyC = atomHop0Key.get(a3)
        hop0KeyD = atomHop0Key.get(a4)
        if hop0KeyA is None or hop0KeyB is None or hop0KeyC is None or hop0KeyD is None:
            raise ValueError(
                f"Improper {impIdx} references missing atom {a1 if hop0KeyA is None else a2 if hop0KeyB is None else a3 if hop0KeyC is None else a4}"
            )

        # Look up improper param
        improperParam = db.lookupImproperParam(hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD)
        if improperParam is None:
            improperNoMatch += 1
            improperNoMatchAtoms.append((a1, a2, a3, a4))
            continue

        # Create improper type if not seen
        coeffs = tuple(improperParam.get("coeffs", []))
        improperTypeId = None
        for itid, existing_coeffs in improperTypeParams.items():
            if len(coeffs) == len(existing_coeffs):
                match = all(abs(c - e) < 0.01 for c, e in zip(coeffs, existing_coeffs))
                if match:
                    improperTypeId = itid
                    break

        if improperTypeId is None:
            improperTypeId = nextImproperTypeId
            nextImproperTypeId += 1
            improperTypeParams[improperTypeId] = coeffs

        improperTypeMap[impIdx] = improperTypeId

    # Write all warnings first
    if bondNoMatch > 0:
        for a1, a2 in bondNoMatchAtoms:
            log.warning(f"  No bond param for atoms {a1}-{a2}")
    if angleNoMatch > 0:
        for a1, a2, a3 in angleNoMatchAtoms:
            log.warning(f"  No angle param for atoms {a1}-{a2}-{a3}")
    if dihedralNoMatch > 0:
        for a1, a2, a3, a4 in dihedralNoMatchAtoms:
            log.warning(f"  No dihedral param for atoms {a1}-{a2}-{a3}-{a4}")
    if improperNoMatch > 0:
        for a1, a2, a3, a4 in improperNoMatchAtoms:
            log.warning(f"  No improper param for atoms {a1}-{a2}-{a3}-{a4}")

    # Write all parameter info lines
    log.info(
        f"Bond parameters: {len(bondTypeMap)} bonds, {len(bondTypeParams)} unique types, {bondNoMatch} no match"
    )
    log.info(
        f"Angle parameters: {len(angleTypeMap)} angles, {len(angleTypeParams)} unique types, {angleNoMatch} no match"
    )
    log.info(
        f"Dihedral parameters: {len(dihedralTypeMap)} dihedrals, {len(dihedralTypeParams)} unique types, {dihedralNoMatch} no match"
    )
    log.info(
        f"Improper parameters: {len(improperTypeMap)} impropers, {len(improperTypeParams)} unique types, {improperNoMatch} no match"
    )

    log.info(
        f"Atom type assignment: {len(atoms)} atoms, {len(atomTypeParams)} unique types, {noMatch} no match"
    )
    log.info(f"  hop3 matches: {hop3Matches}")
    log.info(f"  hop2 matches: {hop2Matches}")
    log.info(f"  hop1 matches: {hop1Matches}")
    log.info(f"  hop0 matches: {hop0Matches}")
    log.info(f"  no match: {noMatch}")

    # Build output data
    # Get box dimensions from coordinates
    allCoords = list(coords.values())
    xvals = [c[0] for c in allCoords]
    yvals = [c[1] for c in allCoords]
    zvals = [c[2] for c in allCoords]
    xlo, xhi = min(xvals) - 5, max(xvals) + 5
    ylo, yhi = min(yvals) - 5, max(yvals) + 5
    zlo, zhi = min(zvals) - 5, max(zvals) + 5

    # Build LammpsData object
    lmpData = LammpsData()
    lmpData.header_comment = "LAMMPS data file Generated by mappingff"
    lmpData.atoms = len(atoms)
    lmpData.bonds = len(bondTypeMap)
    lmpData.angles = len(angleTypeMap)
    lmpData.dihedrals = len(dihedralTypeMap)
    lmpData.impropers = len(improperTypeMap)
    lmpData.atom_types = len(atomTypeParams)
    lmpData.bond_types = len(bondTypeParams)
    lmpData.angle_types = len(angleTypeParams)
    lmpData.dihedral_types = len(dihedralTypeParams)
    lmpData.improper_types = len(improperTypeParams)
    lmpData.xlo = xlo
    lmpData.xhi = xhi
    lmpData.ylo = ylo
    lmpData.yhi = yhi
    lmpData.zlo = zlo
    lmpData.zhi = zhi

    # Masses
    for db_type, params in sorted(atomTypeParams.items()):
        lmpData.masses.append((params["type"], params["mass"]))

    # Pair Coeffs
    for db_type, params in sorted(atomTypeParams.items()):
        lmpData.pair_coeffs.append((params["type"], params["epsilon"], params["sigma"]))

    # Bond Coeffs
    for btid in sorted(bondTypeParams.keys()):
        k, r0 = bondTypeParams[btid]
        lmpData.bond_coeffs.append((btid, k, r0))

    # Angle Coeffs
    for atid in sorted(angleTypeParams.keys()):
        k, theta0 = angleTypeParams[atid]
        lmpData.angle_coeffs.append((atid, k, theta0))

    # Dihedral Coeffs
    for dtid in sorted(dihedralTypeParams.keys()):
        coeffs = dihedralTypeParams[dtid]
        lmpData.dihedral_coeffs.append((dtid,) + coeffs)

    # Improper Coeffs
    for itid in sorted(improperTypeParams.keys()):
        coeffs = improperTypeParams[itid]
        lmpData.improper_coeffs.append((itid,) + coeffs)

    # Atom records
    for atom in atoms:
        atomIdx = atom["idx"]
        lammpsType = atomTypeMap[atomIdx]
        outputType = atomTypeParams[lammpsType]["type"]
        charge = atomTypeParams[lammpsType].get("charge", 0.0)
        x, y, z = coords[atomIdx]
        lmpData.atom_records.append((atomIdx, 1, outputType, charge, x, y, z))

    # Bond records
    for bond in bonds:
        bondIdx = bond["idx"]
        if bondIdx in bondTypeMap:
            bt = bondTypeMap[bondIdx]
            a1 = bond["a1"]
            a2 = bond["a2"]
            lmpData.bond_records.append((bondIdx, bt, a1, a2))

    # Angle records
    for angle in angles:
        angleIdx = angle["idx"]
        if angleIdx in angleTypeMap:
            at = angleTypeMap[angleIdx]
            a1 = angle["a1"]
            a2 = angle["a2"]
            a3 = angle["a3"]
            lmpData.angle_records.append((angleIdx, at, a1, a2, a3))

    # Dihedral records
    for dihedral in dihedrals:
        dihIdx = dihedral["idx"]
        if dihIdx in dihedralTypeMap:
            dt = dihedralTypeMap[dihIdx]
            a1 = dihedral["a1"]
            a2 = dihedral["a2"]
            a3 = dihedral["a3"]
            a4 = dihedral["a4"]
            lmpData.dihedral_records.append((dihIdx, dt, a1, a2, a3, a4))

    # Improper records
    for improper in impropers:
        impIdx = improper["idx"]
        if impIdx in improperTypeMap:
            it = improperTypeMap[impIdx]
            a1 = improper["a1"]
            a2 = improper["a2"]
            a3 = improper["a3"]
            a4 = improper["a4"]
            lmpData.improper_records.append((impIdx, it, a1, a2, a3, a4))

    # Adjust total charge if needed
    before_charge = sum(atom[3] for atom in lmpData.atom_records)
    if total_charge is not None:
        after_step1_charge, after_step2_charge = adjustTotalCharge(
            lmpData, total_charge, db, atomTypeParams
        )
        log.info(
            f"Total charge: {after_step2_charge:.6f} "
            f"(before: {before_charge:.6f}, "
            f"after step1: {after_step1_charge:.6f}, "
            f"after step2: {after_step2_charge:.6f}, "
            f"target: {total_charge})"
        )
    else:
        log.info(f"Total charge: {before_charge:.6f} (no adjustment requested)")

    # Write LAMMPS file
    generateLammps(lmpData, outPath)
    log.info(f"Output written to {outPath}")

    return {
        "atoms": len(atoms),
        "bonds": len(bonds),
        "angles": len(angleTypeMap),
        "dihedrals": len(dihedralTypeMap),
        "impropers": len(improperTypeMap),
        "unique_types": len(atomTypeParams),
        "hop3_matches": hop3Matches,
        "hop2_matches": hop2Matches,
        "hop1_matches": hop1Matches,
        "hop0_matches": hop0Matches,
        "no_match": noMatch,
    }
