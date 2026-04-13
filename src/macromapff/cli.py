"""CLI entry point for MacroMapFF.

Provides three commands:
    - build-db: Build parameter database from sample molecules
    - add-samples: Add new samples to existing database
    - parameterize: Generate parameterized LAMMPS files for target molecules
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

from macromapff.db import MacroMapDB
from macromapff.lmp import parseLammps
from macromapff.mol import MolReader, computeHopKeys
from macromapff.utils import USER_DEFAULT_DB_PATH, setupLogging


def buildDb(samplesDir: Path, dbPath: Path, verbose: bool = False) -> dict:
    """Build a parameter database from sample molecules.

    Iterates over subdirectories in samplesDir, each expected to contain
    a .mol (or .pdb) and a .lmp file with the same base name.

    Args:
        samplesDir: Directory containing sample subdirectories.
        dbPath: Path to output the SQLite database file.
        verbose: If True, print detailed progress.

    Returns:
        Dictionary with build statistics.
    """
    log = logging.getLogger("build-db")

    # Find all sample directories
    sampleDirs = sorted([d for d in samplesDir.iterdir() if d.is_dir()])
    if not sampleDirs:
        log.warning(f"No subdirectories found in {samplesDir}")
        return {"samples_count": 0, "atoms_processed": 0}

    # Initialize database
    db = MacroMapDB(dbPath)
    db.load()

    totalAtoms = 0
    totalBonds = 0
    totalAngles = 0
    totalDihedrals = 0
    totalImpropers = 0

    for sampleDir in sampleDirs:
        segName = sampleDir.name

        # Find mol file
        molFiles = list(sampleDir.glob("*.mol")) + list(sampleDir.glob("*.pdb"))
        lmpFiles = list(sampleDir.glob("*.lmp")) + list(sampleDir.glob("*.lammps.lmp"))

        if not molFiles:
            log.warning(f"No .mol/.pdb file found in {sampleDir}, skipping")
            continue
        if not lmpFiles:
            log.warning(f"No .lmp file found in {sampleDir}, skipping")
            continue

        molPath = molFiles[0]
        lmpPath = lmpFiles[0]

        if verbose:
            log.info(f"Processing {segName}: {molPath.name} + {lmpPath.name}")

        # Parse molecule and LAMMPS data
        molReader = MolReader(molPath)
        lmpData = parseLammps(lmpPath)

        atoms = molReader.getAtoms()
        bonds = molReader.getBonds()
        angles = molReader.getAngles()
        dihedrals = molReader.getDihedrals()
        impropers = molReader.getImpropers()
        coords = molReader.getCoords()

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
            hop2Key, hop1Key, hop0Key = computeHopKeys(molReader, atomIdx)
            hop2Env = molReader.computeHop2Env(atomIdx)
            atomHop0Key[atomIdx] = hop0Key

            # Get LAMMPS type from the sample data
            atomRecord = next(
                (rec for rec in lmpData.atom_records if rec[0] == atomIdx),
                None,
            )
            if atomRecord is None:
                log.warning(f"Atom {atomIdx} not found in LAMMPS data for {segName}")
                continue

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

            # Insert atom type at hop2 level
            # lammps_type is auto-generated to ensure uniqueness per hop2Key
            db.insertAtomType(hop2Key, {
                "element": atom["symbol"],
                "hop1_key": hop1Key,
                "hop0_key": hop0Key,
                "hop2_env": hop2Env,
                "mass": mass[1] if mass else 0.0,
                "sigma": pairCoeff[2] if pairCoeff else 0.0,
                "epsilon": pairCoeff[1] if pairCoeff else 0.0,
                "source": [f"{segName}_{atomIdx}"],
            })

            # Get the auto-generated lammps_type from database
            inserted_info = db.getAtomType(hop2Key)
            assigned_lammps_type = inserted_info["lammps_type"]

            # Insert hop1 keymap (for external validation only)
            db.insertHop1Key(hop1Key, assigned_lammps_type)

            # Insert hop0 keymap (for external validation only)
            db.insertHop0Key(hop0Key, assigned_lammps_type)

        # Process bonds
        # bond_records: [bond_id, bond_type, a1, a2]
        for bondRec in lmpData.bond_records:
            bondType = bondRec[1]
            a1 = bondRec[2]
            a2 = bondRec[3]

            hop0KeyA = atomHop0Key.get(a1)
            hop0KeyB = atomHop0Key.get(a2)
            if hop0KeyA is None or hop0KeyB is None:
                continue

            # Get bond coeffs: [type_id, k, r0]
            bondCoeff = next(
                (rec for rec in lmpData.bond_coeffs if rec[0] == bondType),
                None,
            )
            if bondCoeff is None:
                continue

            # Canonical key (order-independent)
            if hop0KeyA <= hop0KeyB:
                key = (hop0KeyA, hop0KeyB)
            else:
                key = (hop0KeyB, hop0KeyA)

            db.insertBondParam(key, {"k": bondCoeff[1], "r0": bondCoeff[2]})

        # Process angles
        # angle_records: [angle_id, angle_type, a1, a2, a3]
        for angleRec in lmpData.angle_records:
            angleType = angleRec[1]
            a1 = angleRec[2]
            a2 = angleRec[3]
            a3 = angleRec[4]

            hop0KeyA = atomHop0Key.get(a1)
            hop0KeyB = atomHop0Key.get(a2)
            hop0KeyC = atomHop0Key.get(a3)
            if hop0KeyA is None or hop0KeyB is None or hop0KeyC is None:
                continue

            # Get angle coeffs: [type_id, k, theta0]
            angleCoeff = next(
                (rec for rec in lmpData.angle_coeffs if rec[0] == angleType),
                None,
            )
            if angleCoeff is None:
                continue

            # Canonical key (outer atoms swappable)
            if hop0KeyA <= hop0KeyC:
                key = (hop0KeyA, hop0KeyB, hop0KeyC)
            else:
                key = (hop0KeyC, hop0KeyB, hop0KeyA)

            db.insertAngleParam(key, {"k": angleCoeff[1], "theta0": angleCoeff[2]})

        # Process dihedrals
        # dihedral_records: [dih_id, dih_type, a1, a2, a3, a4]
        for dihRec in lmpData.dihedral_records:
            dihType = dihRec[1]
            a1 = dihRec[2]
            a2 = dihRec[3]
            a3 = dihRec[4]
            a4 = dihRec[5]

            hop0KeyA = atomHop0Key.get(a1)
            hop0KeyB = atomHop0Key.get(a2)
            hop0KeyC = atomHop0Key.get(a3)
            hop0KeyD = atomHop0Key.get(a4)
            if hop0KeyA is None or hop0KeyB is None or hop0KeyC is None or hop0KeyD is None:
                continue

            # Get dihedral coeffs: [type_id, k0, k1, k2, k3]
            dihCoeff = next(
                (rec for rec in lmpData.dihedral_coeffs if rec[0] == dihType),
                None,
            )
            if dihCoeff is None:
                continue

            # Canonical key (outer atoms swappable)
            if hop0KeyA <= hop0KeyD:
                key = (hop0KeyA, hop0KeyB, hop0KeyC, hop0KeyD)
            else:
                key = (hop0KeyD, hop0KeyC, hop0KeyB, hop0KeyA)

            db.insertDihedralParam(key, {"coeffs": dihCoeff[1:]})

        # Process impropers
        # improper_records: [imp_id, imp_type, a1, a2, a3, a4]
        for impRec in lmpData.improper_records:
            impType = impRec[1]
            a1 = impRec[2]
            a2 = impRec[3]
            a3 = impRec[4]
            a4 = impRec[5]

            hop0KeyA = atomHop0Key.get(a1)
            hop0KeyB = atomHop0Key.get(a2)
            hop0KeyC = atomHop0Key.get(a3)
            hop0KeyD = atomHop0Key.get(a4)
            if hop0KeyA is None or hop0KeyB is None or hop0KeyC is None or hop0KeyD is None:
                continue

            # Get improper coeffs: [type_id, ...]
            impCoeff = next(
                (rec for rec in lmpData.improper_coeffs if rec[0] == impType),
                None,
            )
            if impCoeff is None:
                continue

            # Canonical key (center fixed, others sorted)
            others = sorted([hop0KeyB, hop0KeyC, hop0KeyD])
            key = (hop0KeyA, others[0], others[1], others[2])

            db.insertImproperParam(key, {"coeffs": impCoeff[1:]})

        if verbose:
            log.info(f"  {segName}: {len(atoms)} atoms, {len(bonds)} bonds, "
                     f"{len(angles)} angles, {len(dihedrals)} dihedrals, {len(impropers)} impropers")

    # Update metadata
    db.setMeta("built_at", datetime.now().isoformat())
    db.setMeta("sample_count", str(len(sampleDirs)))

    # Save database
    db.save()

    log.info(f"Database saved to {dbPath}")
    log.info(f"Samples: {len(sampleDirs)}")
    log.info(f"Atoms/Bonds processed: {totalAtoms} atoms, {totalBonds} bonds, "
             f"{totalAngles} angles, {totalDihedrals} dihedrals, {totalImpropers} impropers")
    log.info(f"Atom type entries (hop2 level): {len(db.atomTypes)}")
    log.info(f"Hop1 keymap entries: {len(db.hop1Keymap)}")
    log.info(f"Hop0 keymap entries: {len(db.hop0Keymap)}")
    log.info(f"  Note: hop0 groups by immediate neighbor signatures only (coarse)")
    log.info(f"  hop1 adds 2nd-order neighbor details (finer)")
    log.info(f"  Therefore hop1 >= hop0 in entry count")
    log.info(f"Bonded parameter entries:")
    log.info(f"  Bonds: {len(db.bondParams)}")
    log.info(f"  Angles: {len(db.angleParams)}")
    log.info(f"  Dihedrals: {len(db.dihedralParams)}")
    log.info(f"  Impropers: {len(db.improperParams)}")

    return {
        "samples_count": len(sampleDirs),
        "atoms_processed": totalAtoms,
        "bonds_processed": totalBonds,
        "angles_processed": totalAngles,
        "dihedrals_processed": totalDihedrals,
        "impropers_processed": totalImpropers,
    }


def parameterize(
    molPath: Path,
    dbPath: Path,
    outPath: Path | None = None,
    verbose: bool = False,
) -> dict:
    """Parameterize a target molecule using the database.

    Args:
        molPath: Path to target molecule .mol or .pdb file.
        dbPath: Path to the SQLite database file.
        outPath: Output LAMMPS file path. If None, uses <mol_file>_param.lmp.
        verbose: If True, print detailed progress.

    Returns:
        Dictionary with parameterization statistics.
    """
    log = logging.getLogger("parameterize")

    # Load database
    db = MacroMapDB(dbPath)
    db.load()
    log.info(f"Database loaded from {dbPath}")
    log.info(f"  Atom types: {len(db.atomTypes)}")
    log.info(f"  Hop1 keymap: {len(db.hop1Keymap)}")
    log.info(f"  Hop0 keymap: {len(db.hop0Keymap)}")
    log.info(f"  Bond params: {len(db.bondParams)}")
    log.info(f"  Angle params: {len(db.angleParams)}")
    log.info(f"  Dihedral params: {len(db.dihedralParams)}")
    log.info(f"  Improper params: {len(db.improperParams)}")

    # Parse target molecule
    molReader = MolReader(molPath)
    atoms = molReader.getAtoms()
    bonds = molReader.getBonds()
    coords = molReader.getCoords()

    log.info(f"Target molecule: {molPath.name}")
    log.info(f"  Atoms: {len(atoms)}")
    log.info(f"  Bonds: {len(bonds)}")

    # Determine output path
    if outPath is None:
        outPath = molPath.with_name(f"{molPath.stem}_param.lmp")

    # Resolve atom types
    from macromapff.fallback import resolveAtomType

    atomTypeMap: dict[int, int] = {}  # atomIdx -> lammpsType (resolved from db)
    atomHop0Key: dict[int, str] = {}  # atomIdx -> hop0Key (resolved from db)
    atomHop2Key: dict[int, str] = {}  # atomIdx -> hop2Key (for typeInfo lookup)
    hop2Matches = 0
    hop1Matches = 0
    hop0Matches = 0
    noMatch = 0

    for atom in atoms:
        atomIdx = atom["idx"]
        hop2Key, hop1Key, hop0Key = computeHopKeys(molReader, atomIdx)
        atomHop0Key[atomIdx] = hop0Key
        atomHop2Key[atomIdx] = hop2Key

        lammpsType, resolvedHop0Key = resolveAtomType(
            hop2Key, hop1Key, hop0Key, atom["symbol"], db
        )
        # Use lammpsType as the type identifier
        atomTypeMap[atomIdx] = lammpsType
        atomHop0Key[atomIdx] = resolvedHop0Key

        # Log atom type assignment
        log.debug(f"  Atom {atomIdx}: element={atom['symbol']}, lammps_type={lammpsType}")

        # Track match statistics
        if hop2Key in db.atomTypes:
            hop2Matches += 1
        elif hop1Key in db.hop1Keymap:
            hop1Matches += 1
        elif hop0Key in db.hop0Keymap:
            hop0Matches += 1
        else:
            noMatch += 1

    # Look up bond parameters
    bondParamMap: dict[int, dict] = {}  # bondIdx -> {k, r0}
    bondTypeMap: dict[int, int] = {}   # bondIdx -> newBondTypeId
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
            bondNoMatch += 1
            bondNoMatchAtoms.append((a1, a2))
            continue

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
            bondTypeParams[bondTypeId] = (paramKey[0], paramKey[1])

        bondParamMap[bondIdx] = bondParam
        bondTypeMap[bondIdx] = bondTypeId

    # Look up angle parameters
    angles = molReader.getAngles()
    angleParamMap: dict[int, dict] = {}
    angleTypeMap: dict[int, int] = {}
    angleTypeParams: dict[int, tuple] = {}
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
            angleNoMatch += 1
            angleNoMatchAtoms.append((a1, a2, a3))
            continue

        # Look up angle param
        angleParam = db.lookupAngleParam(hop0KeyA, hop0KeyB, hop0KeyC)
        if angleParam is None:
            angleNoMatch += 1
            angleNoMatchAtoms.append((a1, a2, a3))
            continue

        # Create angle type if not seen
        paramKey = (angleParam["k"], angleParam["theta0"])
        angleTypeId = None
        for atid, (ak, atheta) in angleTypeParams.items():
            if abs(ak - paramKey[0]) < 0.01 and abs(atheta - paramKey[1]) < 0.1:
                angleTypeId = atid
                break

        if angleTypeId is None:
            angleTypeId = nextAngleTypeId
            nextAngleTypeId += 1
            angleTypeParams[angleTypeId] = (paramKey[0], paramKey[1])

        angleParamMap[angleIdx] = angleParam
        angleTypeMap[angleIdx] = angleTypeId

    # Look up dihedral parameters
    dihedrals = molReader.getDihedrals()
    dihedralParamMap: dict[int, dict] = {}
    dihedralTypeMap: dict[int, int] = {}
    dihedralTypeParams: dict[int, tuple] = {}
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
            dihedralNoMatch += 1
            dihedralNoMatchAtoms.append((a1, a2, a3, a4))
            continue

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

        dihedralParamMap[dihIdx] = dihedralParam
        dihedralTypeMap[dihIdx] = dihedralTypeId

    # Look up improper parameters
    impropers = molReader.getImpropers()
    improperParamMap: dict[int, dict] = {}
    improperTypeMap: dict[int, int] = {}
    improperTypeParams: dict[int, tuple] = {}
    nextImproperTypeId = 1
    improperNoMatch = 0
    improperNoMatchAtoms: list[tuple[int, int, int, int]] = []

    for improper in impropers:
        impIdx = improper["idx"]
        a1 = improper["a1"]
        a2 = improper["a2"]
        a3 = improper["a3"]
        a4 = improper["a4"]
        hop0KeyA = atomHop0Key.get(a1)
        hop0KeyB = atomHop0Key.get(a2)
        hop0KeyC = atomHop0Key.get(a3)
        hop0KeyD = atomHop0Key.get(a4)
        if hop0KeyA is None or hop0KeyB is None or hop0KeyC is None or hop0KeyD is None:
            improperNoMatch += 1
            improperNoMatchAtoms.append((a1, a2, a3, a4))
            continue

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

        improperParamMap[impIdx] = improperParam
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

    # Renumber LAMMPS types to be consecutive starting from 1
    # This is REQUIRED for LAMMPS - type IDs must be consecutive
    uniqueTypes = sorted(set(atomTypeMap.values()))
    typeMapping = {old: new for new, old in enumerate(uniqueTypes, 1)}

    # Write all parameter info lines
    log.info(f"Bond parameters: {len(bondParamMap)} bonds, {len(bondTypeParams)} unique types, {bondNoMatch} no match")
    log.info(f"Angle parameters: {len(angleParamMap)} angles, {len(angleTypeParams)} unique types, {angleNoMatch} no match")
    log.info(f"Dihedral parameters: {len(dihedralParamMap)} dihedrals, {len(dihedralTypeParams)} unique types, {dihedralNoMatch} no match")
    log.info(f"Improper parameters: {len(improperParamMap)} impropers, {len(improperTypeParams)} unique types, {improperNoMatch} no match")

    log.info(f"Atom type assignment: {len(atoms)} atoms, {len(uniqueTypes)} unique types, {noMatch} no match")
    log.info(f"  hop2 exact matches: {hop2Matches}")
    log.info(f"  hop1 fallback: {hop1Matches}")
    log.info(f"  hop0 fallback: {hop0Matches}")
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

    # Build atom type info: output_type -> {element, mass, sigma, epsilon, lammps_type}
    # lammps_type is stored so we can trace back which db entry this came from
    typeInfo: dict[int, dict] = {}
    for atom in atoms:
        atomIdx = atom["idx"]
        lammpsType = atomTypeMap[atomIdx]
        outputType = typeMapping[lammpsType]
        if outputType not in typeInfo:
            element = atom["symbol"]

            # Try to find matching atom type info by element + lammps_type
            info = None
            for key, entry in db.atomTypes.items():
                if entry["element"] == element and entry["lammps_type"] == lammpsType:
                    info = entry
                    break

            if info is not None:
                typeInfo[outputType] = {
                    "element": info["element"],
                    "mass": info["mass"],
                    "sigma": info["sigma"],
                    "epsilon": info["epsilon"],
                    "lammps_type": lammpsType,  # Track original for debugging
                }
            else:
                # No match found - set all parameters to 0 and report error
                log.error(f"  Atom {atomIdx} ({element}): no database match, setting all params to 0")
                typeInfo[outputType] = {
                    "element": element,
                    "mass": 0.0,
                    "sigma": 0.0,
                    "epsilon": 0.0,
                    "lammps_type": lammpsType,
                }

    # Build LammpsData object
    from macromapff.lmp import LammpsData

    lmpData = LammpsData()
    lmpData.header_comment = "LAMMPS data file Generated by MacroMapFF"
    lmpData.atoms = len(atoms)
    lmpData.bonds = len(bondParamMap)
    lmpData.angles = len(angleParamMap)
    lmpData.dihedrals = len(dihedralParamMap)
    lmpData.impropers = len(improperParamMap)
    lmpData.atom_types = len(uniqueTypes)
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
    for newType in sorted(typeInfo.keys()):
        info = typeInfo[newType]
        lmpData.masses.append((newType, info["mass"]))

    # Pair Coeffs
    for newType in sorted(typeInfo.keys()):
        info = typeInfo[newType]
        lmpData.pair_coeffs.append((newType, info["epsilon"], info["sigma"]))

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
        outputType = typeMapping[lammpsType]
        x, y, z = coords[atomIdx]
        lmpData.atom_records.append((atomIdx, 1, outputType, 0.0, x, y, z))

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

    # Write LAMMPS file
    from macromapff.lmp import generateLammps
    generateLammps(lmpData, outPath)
    log.info(f"Output written to {outPath}")

    return {
        "atoms": len(atoms),
        "bonds": len(bonds),
        "angles": len(angleParamMap),
        "dihedrals": len(dihedralParamMap),
        "impropers": len(improperParamMap),
        "unique_types": len(uniqueTypes),
        "hop2_matches": hop2Matches,
        "hop1_matches": hop1Matches,
        "hop0_matches": hop0Matches,
        "no_match": noMatch,
    }


def main() -> None:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="MacroMapFF - Molecular force field parameterization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # build-db command
    build_parser = sub.add_parser("build-db", help="Build parameter database from samples")
    build_parser.add_argument(
        "samples_dir",
        type=Path,
        help="Directory containing sample subdirectories with .mol and .lmp files",
    )
    build_parser.add_argument(
        "--db-dir",
        type=Path,
        default=USER_DEFAULT_DB_PATH.parent,
        help="Output directory for database file (default: ./database)",
    )
    build_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress",
    )

    # add-samples command
    add_parser = sub.add_parser("add-samples", help="Add new samples to existing database")
    add_parser.add_argument(
        "samples_dir",
        type=Path,
        help="Directory containing sample subdirectories",
    )
    add_parser.add_argument(
        "--db-path",
        type=Path,
        default=USER_DEFAULT_DB_PATH,
        help="Path to existing database file (default: ./database/db.pkl)",
    )
    add_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress",
    )

    # parameterize command
    param_parser = sub.add_parser("parameterize", help="Parameterize target molecule")
    param_parser.add_argument(
        "mol_file",
        type=Path,
        help="Path to target molecule .mol or .pdb file",
    )
    param_parser.add_argument(
        "--out",
        type=Path,
        help="Output LAMMPS file path (default: <mol_file>_param.lmp)",
    )
    param_parser.add_argument(
        "--db-path",
        type=Path,
        default=USER_DEFAULT_DB_PATH,
        help="Path to database file (default: ./database/db.pkl)",
    )
    param_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress",
    )

    args = parser.parse_args()

    # Setup logging
    logLevel = logging.DEBUG if getattr(args, "verbose", False) else logging.INFO
    setupLogging(logLevel)

    if args.command == "build-db":
        dbPath = args.db_dir / "db.pkl"
        result = buildDb(args.samples_dir, dbPath, args.verbose)
        print(f"Build complete: {result['samples_count']} samples, {result['atoms_processed']} atoms")

    elif args.command == "add-samples":
        print("add-samples not yet implemented")

    elif args.command == "parameterize":
        result = parameterize(
            args.mol_file,
            args.db_path,
            args.out,
            args.verbose,
        )
        print(f"Parameterize complete: {result['atoms']} atoms, "
              f"{result['bonds']} bonds, "
              f"{result['angles']} angles, "
              f"{result['dihedrals']} dihedrals, "
              f"{result['impropers']} impropers, "
              f"{result['unique_types']} types, "
              f"hop2={result['hop2_matches']}, "
              f"hop1={result['hop1_matches']}, "
              f"hop0={result['hop0_matches']}, "
              f"no_match={result['no_match']}")


if __name__ == "__main__":
    main()
