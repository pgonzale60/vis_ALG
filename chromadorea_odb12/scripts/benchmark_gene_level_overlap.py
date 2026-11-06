#!/usr/bin/env python3
import gzip, os, sys
from collections import defaultdict, Counter

ref_species = [
    ('Caenorhabditis_kamaaina', 'GCA_964211945.1'),
    ('Oscheius_myriophilus', 'GCA_964057215.1'),
    ('Litoditis_marina', 'GCA_964036275.1'),
    ('Auanema_rhodensis', 'GCA_964057225.1'),
    ('Acrobeloides_maximus', 'GCA_964212105.1'),
    ('Acanthocheilonema_viteae', 'GCA_984574255.1'),
    ('Koerneria_pararmatus', 'GCA_964056735.1'),
    ('Cooperia_oncophora', 'GCA_965643665.1')
]

legacy_map = {}
with gzip.open('/lustre/scratch122/tol/teams/blaxter/users/pg17/data/synteny/gene2Nigon_busco20200927.tsv.gz', 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            legacy_map[parts[0]] = parts[1]

nigon12_map = {}
tier12_map = {}
with gzip.open('/lustre/scratch122/tol/teams/blaxter/users/pg17/data/synteny/gene2Nigon_busco_chromadorea_odb12.tsv.gz', 'rt') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            nigon12_map[parts[0]] = parts[1]
            tier12_map[parts[0]] = parts[3]

confusion_matrix = defaultdict(Counter)
matched_pairs = []
overlap_by_species = {}

for sp, acc in ref_species:
    f12 = f'/data/tol/data/tol/nematodes/{sp}/analysis/{acc}/pipelines/busco/busco-sanger-tol-0.3.0-20260715/results/busco/chromadorea_odb12/{acc}.chromadorea_odb12.full_table.tsv.gz'
    f10 = f'/data/tol/data/tol/nematodes/{sp}/analysis/{acc}/pipelines/busco/busco-sanger-tol-0.3.0-20260715/results/busco/nematoda_odb10/{acc}.nematoda_odb10.full_table.tsv.gz'
    
    if not (os.path.exists(f12) and os.path.exists(f10)):
        continue
        
    genes12 = []
    with gzip.open(f12, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if p[1] in ('Complete', 'Duplicated'):
                seq = p[2]
                s = min(int(p[3]), int(p[4]))
                e = max(int(p[3]), int(p[4]))
                genes12.append((seq, s, e, p[0]))
                
    genes10 = []
    with gzip.open(f10, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if p[1] in ('Complete', 'Duplicated'):
                seq = p[2]
                s = min(int(p[3]), int(p[4]))
                e = max(int(p[3]), int(p[4]))
                genes10.append((seq, s, e, p[0]))
                
    scaf10 = defaultdict(list)
    for g in genes10:
        scaf10[g[0]].append(g)
        
    sp_matches = 0
    for g12 in genes12:
        seq12, s12, e12, b12 = g12
        for g10 in scaf10[seq12]:
            _, s10, e10, b10 = g10
            if max(s12, s10) <= min(e12, e10):
                sp_matches += 1
                n12 = nigon12_map.get(b12, 'Unk')
                n10 = legacy_map.get(b10, 'Unassigned')
                t12 = tier12_map.get(b12, 'Unk')
                matched_pairs.append((b12, b10, n12, n10, t12, sp))
                if n10 in ('A','B','C','D','E','N','X') and n12 in ('A','B','C','D','E','N','X'):
                    confusion_matrix[n10][n12] += 1
                break
    overlap_by_species[sp] = (len(genes12), len(genes10), sp_matches)

print("=== Physical Genomic Interval Overlaps by Species ===")
for sp, (n12, n10, n_mat) in overlap_by_species.items():
    print(f"{sp}: {n_mat} overlaps (odb12={n12}, odb10={n10})")

print(f"\nTotal matched pairs across species: {len(matched_pairs)}")

unique_pairs = set((p[0], p[1], p[2], p[3], p[4]) for p in matched_pairs)
print(f"Unique (busco12, busco10) overlapping pairs: {len(unique_pairs)}")

uniq_conf = defaultdict(Counter)
for p in unique_pairs:
    b12, b10, n12, n10, t12 = p
    if n10 in ('A','B','C','D','E','N','X') and n12 in ('A','B','C','D','E','N','X'):
        uniq_conf[n10][n12] += 1

print("\n=== Exact 7x7 Confusion Matrix (Rows = nematoda_odb10 2020, Cols = chromadorea_odb12 2026) ===")
nigons = ['A', 'B', 'C', 'D', 'E', 'N', 'X']
print("odb10\\odb12\t" + "\t".join(nigons) + "\tTotal\tConcordance")
total_diag = 0
total_assigned = 0
for r in nigons:
    row_vals = [uniq_conf[r][c] for c in nigons]
    row_sum = sum(row_vals)
    diag = uniq_conf[r][r]
    total_diag += diag
    total_assigned += row_sum
    pct = (diag / row_sum * 100) if row_sum > 0 else 0
    print(f"{r}\t" + "\t".join(str(v) for v in row_vals) + f"\t{row_sum}\t{pct:.1f}%")

overall_pct = total_diag / total_assigned * 100
print(f"TOTAL Gene-Level Concordance: {total_diag} / {total_assigned} ({overall_pct:.2f}%)")

print("\n=== Legacy 2020 Unassigned Rate for chromadorea_odb12 Tiers ===")
b12_status = defaultdict(list)
for p in unique_pairs:
    b12_status[p[0]].append((p[3], p[4]))

tier_unassigned_counts = Counter()
tier_total_counts = Counter()
for b12, leg_list in b12_status.items():
    t12 = leg_list[0][1]
    is_unassigned = all(item[0] == 'Unassigned' for item in leg_list)
    tier_total_counts[t12] += 1
    if is_unassigned:
        tier_unassigned_counts[t12] += 1

for t in sorted(tier_total_counts.keys()):
    tot = tier_total_counts[t]
    unass = tier_unassigned_counts[t]
    print(f"{t}: {unass} / {tot} unassigned in legacy 2020 ({(unass/tot*100):.1f}%)")
