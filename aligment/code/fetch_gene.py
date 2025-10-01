import os
import sys
import subprocess
import zipfile
from pathlib import Path

#obselete function, only saves meta data on genes, not rRNAs or tRNAs
def fetch_gene(accession: str, target_dir: str):
    # Ensure the target directory exists
    os.makedirs(target_dir, exist_ok=True)

    # Path to the zip file inside the target directory
    outfile = os.path.join(target_dir, f"ncbi_dataset_{accession}.zip")
    outdir = os.path.join(target_dir, f"{accession}_dataset")

    try:
        # Step 1: Download using datasets CLI
        print(f"[INFO] Downloading gene dataset for accession: {accession}")
        subprocess.run(
            ["datasets", "download", "gene", "accession", accession, "--filename", outfile],
            check=True
        )

        # Step 2: Extract the zip into the same target directory
        print(f"[INFO] Extracting {outfile} -> {outdir}/")
        with zipfile.ZipFile(outfile, "r") as zip_ref:
            zip_ref.extractall(outdir)

        # Step 3: Remove the zip
        print(f"[INFO] Cleaning up zip file: {outfile}")
        os.remove(outfile)

        print(f"[DONE] Files are available in: {outdir}/")

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] datasets CLI failed: {e}")
    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: python {sys.argv[0]} <gene_accession> <output_directory>")
        sys.exit(1)

    accession = sys.argv[1]
    output_dir = sys.argv[2]

    fetch_gene(accession, output_dir)

#fetch_gene("KC516842.1", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info")