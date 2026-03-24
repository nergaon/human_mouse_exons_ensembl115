import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    # --- load your tables ---
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    version = "HN6"  
    H_M_table_path = main_dir + f"leafcutter_GSE115736_GSE116177/sum_table_{version}.txt"
    H_M = pd.read_csv(H_M_table_path, sep="\t", index_col=0)
    H_H_table_path = main_dir + f"leafcutter_GSE115736_GSE60424/sum_table_{version}.txt"
    H_H = pd.read_csv(H_H_table_path, sep="\t", index_col=0)
    M_M_table_path = main_dir + f"leafcutter_GSE116177_GSE180020/sum_table_{version}.txt"
    M_M = pd.read_csv(M_M_table_path, sep="\t", index_col=0)
    M_M = M_M.rename(columns={"Cd8T": "CD8T", "BCell": "NveB"})

    dfs = {"H_M": H_M, "H_H": H_H, "M_M": M_M}
    cell_types = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut"]
    # --- extract values ---
    success = {}
    delta02 = {}
    for name, df in dfs.items():
        # reindex ensures missing columns (like Neut) become NaN instead of error
        success[name] = (
        df.loc["Leafcutter success clusters"]
        .reindex(cell_types)
        .astype(float))
        delta02[name] = (
            df.loc["Leafcutter sig and deltapsi>=0.2"]
            .reindex(cell_types)
            .astype(float))
    # --- plotting ---
    x = np.arange(len(cell_types))
    width = 0.25
    fig, ax = plt.subplots(figsize=(12,6))
    for i, (name, _) in enumerate(dfs.items()):
        offset = (i - 1) * width
        # bottom bars (success)
        ax.bar(
            x + offset,
            success[name],
            width,
            label=f"{name} success")    
        # stacked top (delta>=0.2)
        ax.bar(
            x + offset,
            delta02[name],
            width,
            bottom=success[name],
            label=f"{name} ΔPSI≥0.2")
    # --- formatting ---
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=45)
    ax.set_ylabel("Number of clusters")
    ax.set_title("LeafCutter clusters per cell type")
    ax.legend()
    plt.tight_layout()
    output_1 = main_dir + f"leafcutter_comparison_bar_{version}.pdf"
    plt.savefig(output_1, bbox_inches='tight')
    plt.show()

    return

if __name__ == "__main__":
    main()