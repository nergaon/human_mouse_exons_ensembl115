import pandas as pd

BASE = "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026"
files = {
    "A_GSE115736_GSE116177": f"{BASE}/leafcutter_GSE115736_GSE116177/clusters_sum_table_HN6.txt",
    "B_GSE115736_GSE60424": f"{BASE}/leafcutter_GSE115736_GSE60424/clusters_sum_table_HN6.txt",
    "C_GSE116177_GSE180020": f"{BASE}/leafcutter_GSE116177_GSE180020/clusters_sum_table_HN6.txt",
}

for tag, path in files.items():
    df = pd.read_csv(path, sep="\t")
    sub = df[df["genes"].astype(str).str.contains(r"\bUSP16\b", na=False)].copy()
    print(f"\n=== {tag}: USP16 rows = {len(sub)} ===")
    if sub.empty:
        continue

    cols = ["cluster", "h_junction", "genes"]
    for c in ["CD4T_p.adjust", "CD8T_p.adjust", "Cd8T_p.adjust", "Mono_p.adjust", "NK_p.adjust", "NveB_p.adjust", "BCell_p.adjust"]:
        if c in sub.columns:
            cols.append(c)

    print(sub[cols].head(12).to_string(index=False))

# Also check final summary workbook
xl = f"{BASE}/leafcutter_GSE115736_GSE116177/unique_sig_clusters_HN6.xlsx"
summary = pd.read_excel(xl, sheet_name="summary")
hits = summary[summary["genes"].astype(str).str.contains(r"\bUSP16\b", na=False)]
print(f"\n=== Final summary: USP16 rows = {len(hits)} ===")
if not hits.empty:
    show_cols = [c for c in ["cluster", "genes", "CD4T", "CD8T", "NveB", "NK", "Mono", "n_sig"] if c in hits.columns]
    print(hits[show_cols].to_string(index=False))
