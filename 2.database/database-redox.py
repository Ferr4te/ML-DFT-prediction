import pandas as pd
import numpy as np

def get_base(files_name):
    name = files_name
    if name.endswith(".xyz"):
        name = name[:-4]   # remove .xyz

    # Split by last underscore: 1a-00_0n_gn into ["1a-00_0n", "gn"]
    parts = name.rsplit("_", 1)
    if len(parts) == 2 and parts[1] in ("gn", "ox", "rd"):
        return parts[0]
    return name

def build_filename(base, tag):
    return f"{base}_{tag}.xyz"

def main():
    input_path ="dataset_with_HOMO_LUMO_gap.xlsx"
    output_path = "dataset_with_redox_HOMO_LUMO.xlsx"

    HARTREE_TO_KCAL_MOL = 627.51      # 1 hartree = 627.51 kcal/mol
    F_kcal_per_mol_per_V = 23.061     # Faraday constant (kcal/mol/V)
    E_SCE = 4.51                      # SCE reference (V)

    df = pd.read_excel(input_path)
    df['FILES'] = df['FILES'].astype(str)

    # Extract base and tag (gn/ox/rd)
    df['base'] = df['FILES'].apply(get_base)
    df['tag'] = df['FILES'].apply(lambda s: s.replace(".xyz", "").rsplit("_", 1)[-1])

    results = []

    for base in df['base'].unique():
        gn_name = build_filename(base, "gn")
        ox_name = build_filename(base, "ox")
        rd_name = build_filename(base, "rd")

        gn_row = df[df['FILES'] == gn_name]
        ox_row = df[df['FILES'] == ox_name]
        rd_row = df[df['FILES'] == rd_name]

        # if any missing, record NaN
        if gn_row.empty or ox_row.empty or rd_row.empty:
            results.append((base, np.nan))
            continue

        G_gn = float(gn_row.iloc[0]['G (hartree)'])
        G_ox = float(ox_row.iloc[0]['G (hartree)'])
        G_rd = float(rd_row.iloc[0]['G (hartree)'])

        # 1. Convert to kcal/mol
        G_gn_kcal = G_gn * HARTREE_TO_KCAL_MOL
        G_ox_kcal = G_ox * HARTREE_TO_KCAL_MOL
        G_rd_kcal = G_rd * HARTREE_TO_KCAL_MOL

        # 2. Reaction Gibbs free energies (kcal/mol)
        dG_red = G_rd_kcal - G_gn_kcal   # R + e- -> R-
        dG_ox  = G_ox_kcal - G_gn_kcal   # R -> R+ + e-

        # 3. ΔΔrG (kcal/mol)
        delta_delta_G_kcal = dG_red - dG_ox

        # 4. Redox potential (V)
        E_half_V = - (delta_delta_G_kcal / F_kcal_per_mol_per_V) - E_SCE

        results.append((base, E_half_V))

    res_df = pd.DataFrame(results, columns=['base', 'redox_potential_V'])
    df = df.merge(res_df, on='base', how='left')
    df.loc[df['tag'] != "gn", 'redox_potential_V'] = np.nan # Keep groundstate redox value only

    df.to_excel(output_path, index=False)
    print("Saved output to:", output_path)

main()