#!/usr/bin/env python3

import re
import sys
import pandas as pd

def parse_hbond_file(filename):
    parsed_data = []
    # Corrected regex to match input format without explicit bond type
    hbond_pattern = re.compile(
        r"/([AB])\s+(\w+)\s+(\d+)\s+([\w\d']+)\s+/([AB])\s+(\w+)\s+(\d+)\s+([\w\d']+)\s+/([AB])\s+(\w+)\s+(\d+)\s+([\w\d']+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
    )

    with open(filename, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines, 1):
        line = line.strip()
        if not line or '/' not in line or any(s in line for s in ['H-bonds', 'Constraints', 'Models', 'Finding']):
            print(f"Skipping line {i}: {line}")
            continue
        # Normalize multiple spaces or tabs to single space
        line = re.sub(r'\s+', ' ', line)
        match = hbond_pattern.match(line)
        if match:
            donor_chain, donor_res, donor_pos, donor_atom, \
            acceptor_chain, acceptor_res, acceptor_pos, acceptor_atom, \
            h_chain, h_res, h_pos, h_atom, dist, strength = match.groups()
            print(f"Matched line {i}: {line}")
            # Assign aptamer (A) and protein (B) based on chain
            if donor_chain == 'A':
                apt_res, apt_pos, apt_atom = donor_res, donor_pos, donor_atom
                prot_res, prot_pos, prot_atom = acceptor_res, acceptor_pos, acceptor_atom
            else:
                apt_res, apt_pos, apt_atom = acceptor_res, acceptor_pos, acceptor_atom
                prot_res, prot_pos, prot_atom = donor_res, donor_pos, donor_atom

            parsed_data.append({
                "Aptamer Nucleotide": apt_res,
                "Aptamer Position": int(apt_pos),
                "Aptamer Atom": apt_atom,
                "Protein Residue": prot_res,
                "Protein Position": int(prot_pos),
                "Protein Atom": prot_atom,
                "Distance (Å)": float(dist),
                "Bond Strength": float(strength)
            })
        else:
            print(f"Line {i} not matched: {line}")

    df = pd.DataFrame(parsed_data)
    if df.empty:
        print("Warning: No H-bonds parsed. Check input file format or regex pattern.")
    else:
        print(f"Successfully parsed {len(df)} H-bonds.")
    return df

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: ./ChimeraX_hbonds.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    df = parse_hbond_file(input_file)
    output_file = input_file.replace(".txt", "_parsed.csv")
    df.to_csv(output_file, index=False)
    print(f"Data saved to: {output_file}")
