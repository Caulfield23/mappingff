#!/usr/bin/env python3
import csv
import json
from pathlib import Path


def _fmt_box(lo, hi):
    """Format one simulation box axis range for LAMMPS output."""
    return f"{lo:12.6f} {hi:12.6f}"


def write_atom_env_csv(out_dir: Path, module: str, atom_rows):
    """Write one module's atom environment table to CSV."""
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / f"{module}_atom_env.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "module",
                "atom_index",
                "atom_name",
                "opls_type_id",
                "opls_type_name",
                "charge",
                "sigma",
                "epsilon",
                "env_key",
            ],
        )
        writer.writeheader()
        for row in atom_rows:
            writer.writerow(row)

    return csv_path


def write_atom_keytype_map(out_path: Path, atom_records):
    """Write atom-to-key-type mapping for debugging and inspection."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["atom_index", "chosen_key_type", "all_key_types"])
        for rec in atom_records:
            all_key_types = sorted({int(x) for x in rec.get("global_key_ids", [])})
            writer.writerow(
                [
                    int(rec["atom_id"]),
                    int(rec.get("global_key_id", 0)),
                    ";".join(str(x) for x in all_key_types),
                ]
            )


def write_keymap_csv(path: Path, rows, fieldnames):
    """Write generic keymap-like dict rows to a CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_observed_csv(csv_path: Path, observed_list):
    """Write observed multi-atom mapping rows to CSV."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "interaction_kind",
                "n_atoms",
                "lmp_type_tuple",
                "term_count",
                "source_term_types",
                "coeff_param_sets",
            ],
        )
        writer.writeheader()
        for row in observed_list:
            writer.writerow(
                {
                    "interaction_kind": row["interaction_kind"],
                    "n_atoms": row["n_atoms"],
                    "lmp_type_tuple": json.dumps(
                        row["lmp_type_tuple"], ensure_ascii=False, separators=(",", ":")
                    ),
                    "term_count": row["term_count"],
                    "source_term_types": json.dumps(
                        row["source_term_types"],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                    "coeff_param_sets": json.dumps(
                        row["coeff_param_sets"],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                }
            )


def write_master_csv(path: Path, rows):
    """Write merged multi-atom master key-type mapping CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "interaction_kind",
        "key_type_tuple",
        "coeff_param_sets",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_lammps_data(
    out_path: Path,
    mol,
    atom_records,
    atom_type_rows,
    bond_records,
    angle_records,
    dihedral_records,
    improper_records,
    bond_type_rows,
    angle_type_rows,
    dihedral_type_rows,
    improper_type_rows,
    molecule_id: int,
    box_padding: float,
):
    """Serialize assigned parameters and topology back to LAMMPS data format."""
    conf = mol.GetConformer()
    xs, ys, zs = [], [], []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(float(p.x))
        ys.append(float(p.y))
        zs.append(float(p.z))

    xlo, xhi = min(xs) - box_padding, max(xs) + box_padding
    ylo, yhi = min(ys) - box_padding, max(ys) + box_padding
    zlo, zhi = min(zs) - box_padding, max(zs) + box_padding

    atom_by_id = {a["atom_id"]: a for a in atom_records}

    improper_has_zero = any(int(r.get("type_id", 0)) == 0 for r in improper_records)
    if improper_has_zero:
        improper_records_out = []
        for rec in improper_records:
            rec_new = dict(rec)
            rec_new["type_id"] = max(1, int(rec_new.get("type_id", 0)))
            improper_records_out.append(rec_new)
        improper_type_rows_out = improper_type_rows
    else:
        improper_records_out = improper_records
        improper_type_rows_out = improper_type_rows

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("LAMMPS data file generated by parameterize.py\n\n")
        f.write(f"{len(atom_records):8d} atoms\n")
        f.write(f"{len(bond_records):8d} bonds\n")
        f.write(f"{len(angle_records):8d} angles\n")
        f.write(f"{len(dihedral_records):8d} dihedrals\n")
        f.write(f"{len(improper_records_out):8d} impropers\n\n")

        f.write(f"{len(atom_type_rows):8d} atom types\n")
        f.write(f"{len(bond_type_rows):8d} bond types\n")
        f.write(f"{len(angle_type_rows):8d} angle types\n")
        f.write(f"{len(dihedral_type_rows):8d} dihedral types\n")
        f.write(f"{len(improper_type_rows_out):8d} improper types\n\n")

        f.write(f"{_fmt_box(xlo, xhi)} xlo xhi\n")
        f.write(f"{_fmt_box(ylo, yhi)} ylo yhi\n")
        f.write(f"{_fmt_box(zlo, zhi)} zlo zhi\n\n")

        f.write("Masses\n\n")
        for row in atom_type_rows:
            f.write(f"{row['local_type']:8d} {row['mass']:10.6f}\n")
        f.write("\n")

        f.write("Pair Coeffs\n\n")
        for row in atom_type_rows:
            f.write(
                f"{row['local_type']:8d} {row['epsilon']:10.6f} {row['sigma']:10.6f}\n"
            )
        f.write("\n")

        f.write("Bond Coeffs\n\n")
        for tid, coeff in bond_type_rows:
            f.write(f"{tid:8d} " + " ".join(f"{float(x):10.4f}" for x in coeff) + "\n")
        f.write("\n")

        f.write("Angle Coeffs\n\n")
        for tid, coeff in angle_type_rows:
            f.write(f"{tid:8d} " + " ".join(f"{float(x):10.3f}" for x in coeff) + "\n")
        f.write("\n")

        f.write("Dihedral Coeffs\n\n")
        for tid, coeff in dihedral_type_rows:
            f.write(f"{tid:8d} " + " ".join(f"{float(x):10.3f}" for x in coeff) + "\n")
        f.write("\n")

        f.write("Improper Coeffs\n\n")
        for tid, coeff in improper_type_rows_out:
            f.write(
                f"{tid:8d} {float(coeff[0]):10.3f} {int(float(coeff[1])):7d} {int(float(coeff[2])):7d}\n"
            )
        f.write("\n")

        f.write("Atoms\n\n")
        for atom_id in range(1, mol.GetNumAtoms() + 1):
            p = conf.GetAtomPosition(atom_id - 1)
            rec = atom_by_id[atom_id]
            f.write(
                f"{atom_id:8d} {molecule_id:6d} {rec['atom_type']:6d} "
                f"{rec['charge']:10.6f} {float(p.x):10.5f} {float(p.y):10.5f} {float(p.z):10.5f}\n"
            )
        f.write("\n")

        f.write("Bonds\n\n")
        for idx, rec in enumerate(bond_records, start=1):
            i, j = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d}\n")
        f.write("\n")

        f.write("Angles\n\n")
        for idx, rec in enumerate(angle_records, start=1):
            i, j, k = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d} {k:6d}\n")
        f.write("\n")

        f.write("Dihedrals\n\n")
        for idx, rec in enumerate(dihedral_records, start=1):
            i, j, k, l = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d} {k:6d} {l:6d}\n")
        f.write("\n")

        f.write("Impropers\n\n")
        for idx, rec in enumerate(improper_records_out, start=1):
            i, j, k, l = rec["atom_ids"]
            f.write(f"{idx:8d} {rec['type_id']:6d} {i:6d} {j:6d} {k:6d} {l:6d}\n")


class LammpsDataWriter:
    """Object wrapper for writing one LAMMPS data target file."""

    def __init__(self, out_path: Path, molecule_id: int, box_padding: float) -> None:
        """Initialize writer with output path and box writing options."""
        self.out_path = out_path
        self.molecule_id = molecule_id
        self.box_padding = box_padding

    def write(
        self,
        mol,
        atom_records,
        atom_type_rows,
        bond_records,
        angle_records,
        dihedral_records,
        improper_records,
        bond_type_rows,
        angle_type_rows,
        dihedral_type_rows,
        improper_type_rows,
    ):
        """Write one fully parameterized molecule to LAMMPS data file."""
        return write_lammps_data(
            out_path=self.out_path,
            mol=mol,
            atom_records=atom_records,
            atom_type_rows=atom_type_rows,
            bond_records=bond_records,
            angle_records=angle_records,
            dihedral_records=dihedral_records,
            improper_records=improper_records,
            bond_type_rows=bond_type_rows,
            angle_type_rows=angle_type_rows,
            dihedral_type_rows=dihedral_type_rows,
            improper_type_rows=improper_type_rows,
            molecule_id=self.molecule_id,
            box_padding=self.box_padding,
        )