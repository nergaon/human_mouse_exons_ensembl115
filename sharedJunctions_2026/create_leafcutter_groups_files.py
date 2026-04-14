#!/usr/bin/env python3
"""Mirror Leafcutter folders and generate groups_file.txt from AS_clusters_value_HN6 columns."""

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

SOURCE_ROOT = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115")
TARGET_ROOT = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
VALUE_FILE = TARGET_ROOT / "AS_clusters_value_HN6.txt"


def read_samples_from_value_file(value_file: Path) -> list[str]:
    with value_file.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
    if not header:
        raise ValueError(f"Empty or invalid value file: {value_file}")

    # The first column is the junction identifier index.
    return [sample for sample in header[1:] if sample]


def parse_dataset(sample_name: str) -> Optional[str]:
    match = re.search(r"_(GSE\d+)$", sample_name)
    return match.group(1) if match else None


def infer_folder_candidates(sample_name: str) -> list[str]:
    base = re.sub(r"_GSE\d+$", "", sample_name)

    if base.startswith("NveB"):
        return ["NveB", "BCell"]
    if base.startswith(("BCell", "Bcells", "B_fo")):
        return ["BCell", "NveB"]

    if base.startswith(("CD4T", "EffCD4T", "MemCD4T", "NveCD4T", "CD4_", "T_4_")):
        return ["CD4T"]

    if base.startswith(("CD8T", "MemCD8T", "NveCd8T", "CD8_", "T_8_")):
        return ["CD8T", "Cd8T"]

    if base.startswith(("Mono", "InfMono", "Monocytes", "Mo_6C+")):
        return ["Mono"]

    if base.startswith(("Neut", "Neutrophils")):
        return ["Neut"]

    if base.startswith("NK"):
        return ["NK"]

    return []


def get_leafcutter_dirs(source_root: Path) -> list[Path]:
    return sorted(
        d for d in source_root.iterdir()
        if d.is_dir() and d.name.startswith("leafcutter_GSE")
    )


def extract_pair_datasets(leaf_dir_name: str) -> tuple[str, str]:
    match = re.match(r"leafcutter_(GSE\d+)_(GSE\d+)$", leaf_dir_name)
    if not match:
        raise ValueError(f"Unexpected Leafcutter folder name: {leaf_dir_name}")
    return match.group(1), match.group(2)


def ensure_matching_subfolders(source_leaf_dir: Path, target_leaf_dir: Path) -> set[str]:
    group_file_folders: set[str] = set()
    for sub in source_leaf_dir.iterdir():
        if not sub.is_dir():
            continue
        (target_leaf_dir / sub.name).mkdir(parents=True, exist_ok=True)
        if (sub / "groups_file.txt").exists():
            group_file_folders.add(sub.name)
    return group_file_folders


def build_groups_for_leaf_dir(
    samples: list[str],
    dataset_a: str,
    dataset_b: str,
    group_folders: set[str],
) -> dict[str, list[tuple[str, str]]]:
    groups = defaultdict(list)
    seen = set()

    for sample in samples:
        dataset = parse_dataset(sample)
        if dataset not in {dataset_a, dataset_b}:
            continue

        candidates = infer_folder_candidates(sample)
        folder = next((name for name in candidates if name in group_folders), None)
        if folder is None:
            continue

        key = (folder, sample)
        if key in seen:
            continue
        seen.add(key)
        groups[folder].append((sample, dataset))

    return groups


def write_groups_file(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for sample, dataset in rows:
            writer.writerow([sample, dataset])


def main() -> None:
    samples = read_samples_from_value_file(VALUE_FILE)
    leaf_dirs = get_leafcutter_dirs(SOURCE_ROOT)

    if not leaf_dirs:
        raise ValueError(f"No leafcutter_GSE* directories found in {SOURCE_ROOT}")

    for source_leaf_dir in leaf_dirs:
        dataset_a, dataset_b = extract_pair_datasets(source_leaf_dir.name)
        target_leaf_dir = TARGET_ROOT / source_leaf_dir.name
        target_leaf_dir.mkdir(parents=True, exist_ok=True)

        group_folders = ensure_matching_subfolders(source_leaf_dir, target_leaf_dir)
        groups = build_groups_for_leaf_dir(samples, dataset_a, dataset_b, group_folders)

        for folder in sorted(group_folders):
            rows = groups.get(folder, [])
            groups_file_path = target_leaf_dir / folder / "groups_file.txt"
            write_groups_file(groups_file_path, rows)

        print(
            f"Prepared {target_leaf_dir} with {len(group_folders)} group folders "
            f"for {dataset_a} vs {dataset_b}."
        )


if __name__ == "__main__":
    main()
