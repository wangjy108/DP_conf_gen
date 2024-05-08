import os
import logging
from pathlib import Path
import argparse
import shutil


def main():
    parser = argparse.ArgumentParser(description='topo I binding pose from NLdock')
    parser.add_argument('--input_csv', type=str, required=True, 
                        help="input csv with molecule smiles, \
                            column with smiles should have header name [SMILES],\
                            column with molecule name should have header name [NAME]")
    parser.add_argument('--energy_window', type=float, default=5.0,
                        help="energy window applied for conformation sampling, unit in kcal/mol, default is 5.0")
    parser.add_argument('--rmsd_cutoff', type=float, default=1.5,
                        help="rmsd cutoff applied for conformation sampling, unit in Angstrom, default is 1.5")
    parser.add_argument('--n_pose_per_conformer', type=int, default=3,
                        help="saved docking pose for each conformer, default is 3")
    
    args = parser.parse_args()

    from pseudo import confGen_topo

    sdf_content, error_msg = confGen_topo(input_csv=args.input_csv,
                                          energy_window=args.energy_window,
                                          rmsd_cutoff=args.rmsd_cutoff,
                                          n_pose_per_conformer=args.n_pose_per_conformer).run()
    
    output_dir = Path("./")

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "output.sdf").open("w+") as f:
        for each in sdf_content:
            f.write(each)

    if error_msg:
        with (output_dir / "ERROR.smi").open("w+") as ff:
            for ee in error_msg:
                ff.write(ee)
    
    logging.info("Check ./outputs for final results")

    os.remove("gen.py")
    os.remove("pseudo.py")
    shutil.rmtree("__pycache__")


if __name__ == '__main__':
    _ = shutil.copy2("/opt/scripts/gen.py", "./")
    _ = shutil.copy2("/opt/scripts/pseudo.py", "./")

    main()
    