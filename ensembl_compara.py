import requests
import pandas as pd

ENSEMBL_REST = "https://rest.ensembl.org"

def get_exons(transcript_id):
    """Return exon coordinates for any Ensembl transcript (human, mouse, etc.)."""
    # Detect species from prefix
    if transcript_id.startswith("ENSMUST"):
        species = "mouse"
    elif transcript_id.startswith("ENST"):
        species = "human"
    else:
        species = None

    url = f"{ENSEMBL_REST}/lookup/id/{transcript_id}"
    params = {"expand": "1"}
    if species and species != "human":
        params["species"] = species  # required for mouse and others
    headers = {"Content-Type": "application/json"}

    r = requests.get(url, params=params, headers=headers)
    r.raise_for_status()
    data = r.json()

    if "Exon" not in data:
        raise ValueError(f"No exons found for {transcript_id}")

    exons = [
        {
            "exon_id": e["id"],
            "seq_region": e["seq_region_name"],
            "start": e["start"],
            "end": e["end"],
            "strand": e["strand"],
        }
        for e in data["Exon"]
    ]
    return pd.DataFrame(exons)

def map_region_to_mouse(chrom, start, end):
    """Use Ensembl Compara to map a human region to mouse coordinates."""
    url = f"{ENSEMBL_REST}/map/human/GRCh38/{chrom}:{start}..{end}/mouse"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers)
    r.raise_for_status()
    mappings = r.json().get("mappings", [])
    mapped = [
        {
            "mouse_seq_region": m["mapped"]["seq_region_name"],
            "mouse_start": m["mapped"]["start"],
            "mouse_end": m["mapped"]["end"],
            "mouse_strand": m["mapped"]["strand"],
        }
        for m in mappings
    ]
    return mapped


def find_overlaps(human_df, mouse_df, overlap_threshold=0.5):
    """Identify mouse exons overlapping mapped human exons."""
    results = []
    for _, row in human_df.iterrows():
        #print(row)
        mappings = map_region_to_mouse(row.seq_region, row.start, row.end)
        for m in mappings:
            overlaps = mouse_df[
                (mouse_df.seq_region == m["mouse_seq_region"]) &
                (mouse_df.end >= m["mouse_start"]) &
                (mouse_df.start <= m["mouse_end"])
            ]
            for _, mo in overlaps.iterrows():
                overlap_len = min(m["mouse_end"], mo.end) - max(m["mouse_start"], mo.start)
                exon_len = mo.end - mo.start
                frac = max(0, overlap_len) / exon_len
                if frac >= overlap_threshold:
                    results.append({
                        "human_exon_id": row.exon_id,
                        "human_coord": f"{row.seq_region}:{row.start}-{row.end}",
                        "mouse_exon_id": mo.exon_id,
                        "mouse_coord": f"{mo.seq_region}:{mo.start}-{mo.end}",
                        "overlap_fraction": round(frac, 2)
                    })
    return pd.DataFrame(results)


# ==== Example run ====
human_transcript = "ENST00000436547"    # human 
mouse_transcript = "ENSMUST00000106226" # mouse 

print("Fetching exons...")
human_exons = get_exons(human_transcript)
mouse_exons = get_exons(mouse_transcript)

print(f"Human exons: {len(human_exons)}, Mouse exons: {len(mouse_exons)}")

print("Mapping and comparing...")
orthologs = find_overlaps(human_exons, mouse_exons)

print(f"Found {len(orthologs)} probable exon orthologs.")
print(orthologs)