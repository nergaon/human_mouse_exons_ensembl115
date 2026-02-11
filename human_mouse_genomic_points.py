import requests
import pandas as pd

ENSEMBL_REST = "https://rest.ensembl.org"

def map_human_to_mouse(chrom, start, end, strand=1):
    """
    Map a human genomic region to orthologous regions in the mouse genome
    using Ensembl Compara REST API.

    Parameters
    ----------
    chrom : str
        Chromosome name (e.g., "X", "7", "1")
    start : int
        Start coordinate (1-based, inclusive)
    end : int
        End coordinate (1-based, inclusive)
    strand : int
        1 for forward, -1 for reverse (default = 1)

    Returns
    -------
    pd.DataFrame
        DataFrame of orthologous regions in mouse with coordinates and confidence score.
    """
    url = f"{ENSEMBL_REST}/map/human/GRCh38/{chrom}:{start}..{end}:{strand}/mouse"
    headers = {"Content-Type": "application/json"}

    r = requests.get(url, headers=headers)
    r.raise_for_status()
    data = r.json()

    mappings = data.get("mappings", [])
    if not mappings:
        print("⚠️ No orthologous region found for this input.")
        return pd.DataFrame()

    results = []
    for m in mappings:
        mapped = m["mapped"]
        results.append({
            "human_coord": f"{chrom}:{start}-{end}:{strand}",
            "mouse_chr": mapped["seq_region_name"],
            "mouse_start": mapped["start"],
            "mouse_end": mapped["end"],
            "mouse_strand": mapped["strand"],
            "alignment_length": mapped["end"] - mapped["start"] + 1
        })

    return pd.DataFrame(results)

chrom = "10"
start = 119573465
end = 119596964
strand = -1
map_human_to_mouse(chrom, start, end, strand)