import os
from Bio import SeqIO

#import symbol key and symbol formatting function
from symbol_dictionary import symbol_correction
from symbol_dictionary import symbol_key as symbol_key

#read data file in gbk format
def read_gbk(accession: str, data_folder: str):
    """
    Reads a GenBank (.gbk) file and extracts coding genes into a list of dictionaries.
    - Uses 'standard_name' if available; otherwise falls back to 'gene', then 'product'.
    - Excludes any genes whose symbol contains 'trna'.
    - Works for single-gene or multi-gene records.
    
    Returns a list of dicts with keys: symbol, start, end, orientation, species.
    """
    gbk_path = os.path.join(data_folder, f"{accession}.gbk")
    if not os.path.exists(gbk_path):
        raise FileNotFoundError(f"No .gbk file found at: {gbk_path}")

    gene_list = []

    with open(gbk_path, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            species = record.annotations.get("organism", "Unknown")
            for feature in record.features:
                if feature.type in ("CDS", "rRNA", "misc_feature"):  # only coding sequences or rRNAs
                    qualifiers = feature.qualifiers

                    # Try multiple keys for gene symbol
                    symbol = qualifiers.get("standard_name", [None])[0]
                    if not symbol:
                        symbol = qualifiers.get("gene", [None])[0]
                    if not symbol:
                        symbol = qualifiers.get("product", [None])[0]
                    # Skip if no symbol or it's a tRNA gene
                    if not symbol or "trna" in symbol.lower():
                        continue
                    
                    #make sure symbol format is standardized
                    symbol_correct = symbol_correction(symbol, symbol_key)
                    gene_info = {
                        "symbol": symbol_correct,
                        "begin": int(feature.location.start) + 1,  # convert 0-based to 1-based
                        "end": int(feature.location.end),
                        "orientation": "+" if feature.location.strand == 1 else "-",
                        "species": species.strip().replace(" ", "-"),
                        "accession": accession
                    }
                    gene_list.append(gene_info)
    #for gene in gene_list: #check necessary genes are present in returned data
    #    print(gene["symbol"])

    return gene_list

#print(read_gbk("OZ208999.1", "fasta_info"))
#print(read_gbk("KC516842.1", "fasta_info"))
#print(read_gbk("NC_010298", "fasta_info"))