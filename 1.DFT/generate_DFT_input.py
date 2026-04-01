import os

def getMolecuelsInfo(filename):
    chargeDic = {"0n":0, "1a":-1, "2a":-2, "1c":1, "2c":2 }
    base = os.path.basename(filename)
    base = base.replace(".xyz", "")
    parts = base.split("_")
    left = parts[0]              # "family-idx"
    family, idx = left.split("-")
    charge = parts[1]            # "0n", "1a", etc.
    state = parts[2] if len(parts) > 2 else "unknown"
    return family, idx, chargeDic[charge], state, base

def main():
    directory_path = "redox_pbe0"  
    for filename in os.listdir(directory_path):
        if filename.endswith("gn.xyz"):
            file_path = os.path.join(directory_path, filename)
            if os.path.isfile(file_path):
                try:
                    family, idx, charge, state, base = getMolecuelsInfo(filename)
                    output_name = filename.replace(".xyz", ".inp")
                    with open(output_name, 'w', encoding='utf-8') as wf:
                        wf.write(f"!PBE0 def2-tzvpp def2/JK RIJK\n%pal nprocs 8 end\n%maxcore 1000\n* xyz {charge} 1\n")

                        with open(file_path, 'r', encoding='utf-8') as f:
                            content = f.read()
                            lines = content.split("\n")
                            for i in lines[2:]:  # Skip first two lines to get coordinates
                                if i.strip():    # Skip empty lines
                                    wf.write(i + "\n")

                        wf.write("*\n")  

                except Exception as e:
                    print(e)

main()
