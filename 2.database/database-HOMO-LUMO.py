import os
import pandas as pd
import re

def parse_homo_lumo_from_result(path_result_out):
    in_table = False
    orbitals = []  # list of (NO, OCC, E_h)

    with open(path_result_out, "r") as f:
        for line in f:
            # Detect the start of the ORBITAL ENERGIES block
            if line.strip().startswith("ORBITAL ENERGIES"):
                in_table = False  
                continue

            # The header line with NO OCC E(Eh) E(eV)
            if re.search(r"\bNO\s+OCC\s+E\(Eh\)", line):
                in_table = True
                continue

            # Stop if we reached blank line 
            if in_table:
                # End condition: often a line starting with '*Only' or a blank line
                if line.strip().startswith("*Only") or not line.strip():
                    break

                # Expected table row format, e.g.:
                #   63   2.0000      -0.197459        -5.3731
                parts = line.split()
                if len(parts) < 3:
                    continue  

                try:
                    no = int(parts[0])
                    occ = float(parts[1])
                    e_h = float(parts[2])
                    orbitals.append((no, occ, e_h))
                except ValueError:
                    continue

    if not orbitals:
        raise ValueError(f"{path_result_out}")

    # Identify HOMO: last orbital with OCC > 0
    occupied = [orb for orb in orbitals if orb[1] > 0.0]
    virtual  = [orb for orb in orbitals if orb[1] == 0.0]

    if not occupied or not virtual:
        raise ValueError(f"{path_result_out}")

    homo_no, homo_occ, homo_e = occupied[-1]  # last occupied
    lumo_no, lumo_occ, lumo_e = virtual[0]    # first virtual

    gap = lumo_e - homo_e  # Hartree

    return homo_e, lumo_e, gap

def main():
    excel_path = "dataset.xlsx" # Same directory
    results_dir = "results"     # Same directory


    df = pd.read_excel(excel_path, sheet_name="Sheet1")
    # Prepare new columns
    df["HOMO_Eh"] = None
    df["LUMO_Eh"] = None
    df["HOMO_LUMO_gap_Eh"] = None
    df["gap_parse_ok"] = False  

    # Loop over rows and write each corresponding result.out
    for idx, row in df.iterrows():
        xyz_name = row["FILES"]             # e.g. "1a-00_0n_gn.xyz"
        base, _ = os.path.splitext(xyz_name)
        out_name = base + ".out"           # e.g. "1a-00_0n_gn.out"
        out_path = os.path.join(results_dir, out_name)

        if not os.path.isfile(out_path):
            print(f"[WARN] Missing result file for {xyz_name}: {out_path}")
            continue

        try:
            homo_e, lumo_e, gap_e = parse_homo_lumo_from_result(out_path)
            df.at[idx, "HOMO_Eh"] = homo_e
            df.at[idx, "LUMO_Eh"] = lumo_e
            df.at[idx, "HOMO_LUMO_gap_Eh"] = gap_e
            df.at[idx, "gap_parse_ok"] = True
        except Exception as e:
            print(f"[ERROR] Failed to parse {out_path}: {e}")

    output_excel = "dataset_with_HOMO_LUMO_gap.xlsx"
    df.to_excel(output_excel, sheet_name="Sheet1", index=False)

    print(f"Done. Saved: {output_excel}")

main()