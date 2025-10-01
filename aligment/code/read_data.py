import os
import json

#obselete function as JSON file does not include data on rRNAs
def read_data(accession: str, data_folder: str):
    """
    Reads a JSON dataset for a given accession and extracts gene information.
    Returns a list of dictionaries with the same structure as read_gbk():
    - symbol
    - start
    - end
    - orientation (+/-)
    - species (taxname)
    Only includes genes with a standard_name (symbol) and skips tRNA genes.
    """
    # path to the JSONL file inside the dataset folder
    json_path = os.path.join(data_folder, f"{accession}_dataset", "ncbi_dataset", "data", "data_report.jsonl")

    if not os.path.exists(json_path):
        raise FileNotFoundError(f"No data_report.jsonl found at: {json_path}")

    gene_list = []

    with open(json_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)

            # 'obj' can contain multiple genes; check if it's a gene entry
            # Some entries may not have 'symbol' or 'annotations'
            symbol = obj.get("symbol")
            if not symbol or "trna" in symbol.lower():
                continue

            # Get genomic location
            start, end, orientation = None, None, "+"
            annotations = obj.get("annotations", [])
            if annotations:
                # Take the first annotation with a genomic location
                for ann in annotations:
                    locations = ann.get("genomicLocations", [])
                    if locations:
                        loc = locations[0].get("genomicRange", {})
                        start = int(loc.get("begin", 0))
                        end = int(loc.get("end", 0))
                        orientation = loc.get("orientation", "+")
                        break  # only use first location

            species = obj.get("taxname", "Unknown")

            gene_list.append({
                "symbol": symbol,
                "start": start,
                "end": end,
                "orientation": orientation,
                "species": species,
                "accession": accession
            })

    return gene_list

print(read_data("NC_024563.1", "fasta_info"))#[0]["symbol"])