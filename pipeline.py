#!/usr/bin/env python3

import os
import sys
import json
import shutil
import subprocess
import pandas as pd
import re
import glob

# AA mapping
aa_dict = {k: v for k, v in {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASH': 'D', 'CYS': 'C', 'CYX': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLH': 'E', 'GLY': 'G', 'HIS': 'H', 'HIE': 'H', 'HID': 'H', 'HIP': 'H',
    'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}.items()}
solvent_resnames = {"HOH", "WAT", "TIP3"}

def run_command(command, stdout_file=None, cwd=None):
    print(f"Running: {' '.join(command)}")
    try:
        if stdout_file:
            with open(stdout_file, 'w') as f_out:
                subprocess.run(command, check=True, text=True, stdout=f_out, stderr=subprocess.STDOUT, cwd=cwd)
        else:
            subprocess.run(command, check=True, text=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with code {e.returncode}")
        sys.exit(1)

def replace_aa_codes(s):
    substrings = s.split(':')
    return ':'.join(aa_dict.get(substr.upper(), substr) for substr in substrings)

def parse_selection_range(selection_str):
    match = re.search(r"residue\s+(\d+)\s+to\s+(\d+)", selection_str)
    if match:
        start, end = map(int, match.groups())
        range_set = set(map(str, range(start, end+1)))
        range_str = f":{start}-{end}"
        return range_set, range_str
    return set(), ":0-0"

def process_contact_freq(tsv_file, threshold=0.1):
    df = pd.read_csv(tsv_file, names=["atom_1", "atom_2", "freq"], sep='\s+', comment="#")
    df = df[df['freq'] >= threshold]
    df["res_1"] = [":".join(x.split(":")[1:3]) for x in df.atom_1]
    df["res_2"] = [":".join(x.split(":")[1:3]) for x in df.atom_2]
    df['res_1'] = df['res_1'].apply(replace_aa_codes)
    df['res_2'] = df['res_2'].apply(replace_aa_codes)
    code = os.path.splitext(os.path.basename(tsv_file))[0]
    df = df.assign(C=df['res_1'] + '-' + df['res_2'])
    df = df.drop(['atom_1', 'atom_2'], axis=1)
    df = df.rename(columns={'freq': code})
    return df, code

def write_bfactors_to_pdb(input_pdb, output_pdb, bfactors):
    shutil.copy(input_pdb, output_pdb)
    with open(output_pdb, 'r') as f:
        lines = f.readlines()
    with open(output_pdb, 'w') as f:
        for line in lines:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname in solvent_resnames:
                    f.write(line)
                    continue
                resnum = line[22:26].strip()
                if resnum in bfactors:
                    bfac = bfactors[resnum]
                    line = line[:60] + f"{bfac:6.2f}" + line[66:]
            f.write(line)

def filter_high_bfactor_residues(bfactors, threshold, allowed_range):
    return {resnum: b for resnum, b in bfactors.items()
            if float(b) > threshold and resnum in allowed_range}

def write_cpptraj_script(residues1, target_range2_str, residues2, target_range1_str, top, traj, output_script):
    with open(output_script, 'w') as f:
        f.write(f"parm ../../{top}\n")
        f.write(f"trajin ../../{traj}\n")
        f.write("go\n")

        for i, resnum in enumerate(residues1, 1):
            f.write(f"lie res{i} :{resnum} {target_range2_str} out int_{resnum}.dat\n")

        offset = len(residues1) + 1
        for j, resnum in enumerate(residues2, offset):
            f.write(f"lie res{j} :{resnum} {target_range1_str} out int_{resnum}.dat\n")

        f.write("go\nquit\n")

def parse_energy_data(energy_dir):
    ele_data, vdw_data = [], []
    for filepath in glob.glob(os.path.join(energy_dir, "int*.dat")):
        try:
            res_id = os.path.basename(filepath).split("_")[-1].replace(".dat", "")
            df = pd.read_csv(filepath, sep='\s+', names=['FRAME', 'ELE', 'VDW'], comment='#')
            ele_mean, ele_std = df['ELE'].mean(), df['ELE'].std()
            vdw_mean, vdw_std = df['VDW'].mean(), df['VDW'].std()
            ele_data.append([res_id, f"{ele_mean:.2f} ± {ele_std:.2f}"])
            vdw_data.append([res_id, f"{vdw_mean:.2f} ± {vdw_std:.2f}"])
        except Exception as e:
            print(f"Error processing {filepath}: {e}")

    pd.DataFrame(ele_data, columns=["Residue", "ELE (mean ± std)"]).to_csv(os.path.join(energy_dir, "ele_summary.csv"), index=False)
    pd.DataFrame(vdw_data, columns=["Residue", "VDW (mean ± std)"]).to_csv(os.path.join(energy_dir, "vdw_summary.csv"), index=False)
    print(f"→ Energy summaries written to {energy_dir}")

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} config.json")
        sys.exit(1)

    config_path = sys.argv[1]
    with open(config_path) as f:
        config = json.load(f)

    topology = config["topology_file"]
    trajectory = config["trajectory_file"]
    outdir = config.get("output_dir", ".")
    os.makedirs(outdir, exist_ok=True)

    # Step 1: Dynamic Contacts
    if "dynamic_contacts" in config:
        print("\n--- [1/7] Calculating dynamic contacts ---")
        dc = config["dynamic_contacts"]
        tsv_out = os.path.join(outdir, dc["output_tsv"])
        cmd = [
            "get_dynamic_contacts.py",
            "--topology", topology,
            "--trajectory", trajectory,
            "--itypes", *dc["itypes"].split(),
            "--sele", dc["selection1"],
            "--sele2", dc["selection2"],
            "--output", tsv_out
        ]
        run_command(cmd)

    # Step 2: Contact Frequencies
    print("\n--- [2/7] Calculating contact frequencies ---")
    freq_tsv = os.path.join(outdir, config["frequency_analysis"]["output_freq_tsv"])
    tsv_in = os.path.join(outdir, config["dynamic_contacts"]["output_tsv"])
    run_command(["get_contact_frequencies.py", "--input_file", tsv_in, "--output_file", freq_tsv])

    # Step 3: Generate PDB
    print("\n--- [3/7] Generating reference PDB ---")
    pdb_out = os.path.join(outdir, config["pdb_generation"]["output_pdb"])
    run_command(["ambpdb", "-p", topology, "-c", trajectory], stdout_file=pdb_out)

    # Step 4: Map Contact Frequencies to PDB
    print("\n--- [4/7] Mapping contact frequencies to PDB ---")
    df, code = process_contact_freq(freq_tsv)
    bfactors_raw = {}
    for _, row in df.iterrows():
        for res in (row["res_1"], row["res_2"]):
            resnum = res.split(":")[1]
            bfactors_raw[resnum] = bfactors_raw.get(resnum, 0.0) + row[code]
    mapped_pdb = os.path.join(outdir, f"{code}_bfac.pdb")
    write_bfactors_to_pdb(pdb_out, mapped_pdb, bfactors_raw)

    # Step 5: Filter high B-factor residues
    print("\n--- [5/7] Filtering high B-factor residues ---")
    bconfig = config["bfactor_filtering"]
    threshold = bconfig["bfactor_threshold"]

    sel1 = config["dynamic_contacts"]["selection1"]
    sel2 = config["dynamic_contacts"]["selection2"]
    range1_set, range1_str = parse_selection_range(sel1)
    range2_set, range2_str = parse_selection_range(sel2)

    residues1 = filter_high_bfactor_residues(bfactors_raw, threshold, range1_set)
    residues2 = filter_high_bfactor_residues(bfactors_raw, threshold, range2_set)

    print(f"→ Range1 high-B: {list(residues1.keys())}")
    print(f"→ Range2 high-B: {list(residues2.keys())}")

    # Step 6: Write LIE cpptraj script
    print("\n--- [6/7] Writing LIE cpptraj script ---")
    energy_dir = os.path.join(outdir, "energies")
    os.makedirs(energy_dir, exist_ok=True)
    script_path = os.path.join(energy_dir, bconfig["output_script"])
    write_cpptraj_script(
        list(residues1.keys()), range2_str,
        list(residues2.keys()), range1_str,
        topology, trajectory, script_path
    )
    print(f"→ LIE script written to {script_path}")

    # Step 7: Optional energy calculation
    if config.get("calculate_energies", False):
        print("\n--- [7/7] Running cpptraj & summarizing energies ---")
        #run_command(["cpptraj", script_path], cwd=energy_dir)
        run_command(["cpptraj", os.path.basename(script_path)], cwd=energy_dir)
        parse_energy_data(energy_dir)

    print("\nContact analysis complete!")

if __name__ == "__main__":
    main()
