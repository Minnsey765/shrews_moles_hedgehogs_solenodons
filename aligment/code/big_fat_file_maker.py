#function to concatenate all sequences into a single sequence to be aligned at once

import os
import sys
from overlord import overlord_function
from matcher import find_best_key
import pandas as pd
from Bio import SeqIO
from pathlib import Path

def big_fat_file_maker(raw_fol: str, csv_path: str, output: str):

    #make list of all species names
    csv = pd.read_csv(csv_path)
    species_names = csv["Species"].tolist()
    #make list of all gene names (ignore species column)
    gene_names = csv.columns[1:].tolist()

    #convert list to dictionary where entries in list are dictionary keys
    my_dict = {item: "" for item in species_names}
    folder = Path(raw_fol)
    for name in gene_names:
        matching_files = list(folder.glob(f"*{name}*.nex"))
        
        if not matching_files:
            print(f"No NEXUS file found for {name}")
            continue

        for file_path in matching_files:
            print(f"Reading {file_path.name}")
            for record in SeqIO.parse(file_path, "nexus"):
                #ignore extra genes for a single species
                #print(len(record.seq))
                if record.id.endswith(".copy"):
                    continue
                    #print(f">{record.id}")

                #use matching function to locate correct key
                match = find_best_key(record.id, my_dict)
                my_dict[match] += str(record.seq)
                #print(len(my_dict[match]))
    
    os.makedirs(output, exist_ok=True)
    return(my_dict)

    


(big_fat_file_maker("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/aligned_fastas", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/newGenBank.csv", "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data"))