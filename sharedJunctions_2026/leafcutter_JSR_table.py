import argparse
import os
import sys
import gzip
from collections import defaultdict

def parse_junction_file(junction_file, max_intron_len, chrom_lst, nochromcheck, ignore_strand=False):
    by_chrom = defaultdict(lambda: defaultdict(int))
    if junction_file.endswith('.gz'):
        opener = gzip.open(junction_file, 'rt')
    else:
        opener = open(junction_file, 'r')

    with opener as f:
        for ln in f:
            lnsplit = ln.split()
            if len(lnsplit) < 12:
                sys.stderr.write("Error in %s\n" % junction_file)
                continue
            chrom, A, B, dot, counts, strand, rA, rb, rgb, blockCount, blockSize, blockStarts = lnsplit
            if int(blockCount) > 2:
                continue
            if not ignore_strand and strand == "?":
                continue
            if not nochromcheck and (chrom not in chrom_lst):
                continue
            Aoff, Boff = blockSize.split(",")
            A, B = int(A) + int(Aoff), int(B) - int(Boff) + 1
            if B - A > int(max_intron_len):
                continue
            if ignore_strand:
                strand = "."
            try:
                by_chrom[(chrom, strand)][(A, B)] += int(counts)
            except:
                by_chrom[(chrom, strand)] = {(A, B): int(counts)}
    return by_chrom

def main():
    parser = argparse.ArgumentParser(description="Generate junction read count table without clustering.")
    parser.add_argument('-j', '--junction_list', required=True, help="File containing list of junction files (one per line).")
    parser.add_argument('-o', '--output_prefix', required=True, help="Output prefix for the TSV file.")
    parser.add_argument('-l', '--max_intron_len', required=True, help="Maximum intron length.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output.")
    parser.add_argument('--nochromcheck', action='store_true', help="Skip chromosome check.")
    parser.add_argument('--ignore_strand', dest='ignore_strand', action='store_true', default=True,
                        help="Ignore strand information; merge counts across strands (default).")
    parser.add_argument('--keep_strand', dest='ignore_strand', action='store_false',
                        help="Keep strand information in junction IDs (strand-aware output).")
    args = parser.parse_args()

    chrom_lst = ["chr%d" % x for x in range(1, 23)] + ['chrX', 'chrY'] + ["%d" % x for x in range(1, 23)] + ['X', 'Y']
    sample_junction_counts = defaultdict(lambda: defaultdict(int))
    sample_names = []

    with open(args.junction_list, 'r') as f:
        flist = f.readlines()

    for k, libl in enumerate(flist, 1):
        lib = libl.strip()
        if not os.path.isfile(lib):
            raise ValueError("File %s does not exist." % lib)
        if args.verbose:
            sys.stderr.write("scanning %d/%d: %s ...\n" % (k, len(flist), lib))

        sample_name = os.path.basename(lib).split('.')[0]
        sample_names.append(sample_name)
        junction_counts = parse_junction_file(lib, args.max_intron_len, chrom_lst, args.nochromcheck, args.ignore_strand)
        for (chrom, strand), junctions in junction_counts.items():
            for (A, B), count in junctions.items():
                junction_id = f"{chrom}:{A}:{B}:{strand}"
                sample_junction_counts[junction_id][sample_name] = count

    output_file = f"{args.output_prefix}_junction_counts_no_strand.tsv"
    with open(output_file, 'w') as f:
        header = ["junction_id"] + sample_names
        f.write("\t".join(header) + "\n")
        for junction_id, counts in sample_junction_counts.items():
            row = [junction_id] + [str(counts[sample]) for sample in sample_names]
            f.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
