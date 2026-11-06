import os, glob, json

base = '/data/tol/data/tol/nematodes'
species_dirs = sorted([d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d))])

taxonomy_lookup = {
    'Caenorhabditis': ('Rhabditida', 'Rhabditina', 'Rhabditidae'),
    'Oscheius': ('Rhabditida', 'Rhabditina', 'Rhabditidae'),
    'Auanema': ('Rhabditida', 'Rhabditina', 'Rhabditidae'),
    'Litoditis': ('Rhabditida', 'Rhabditina', 'Rhabditidae'),
    'Diploscapter': ('Rhabditida', 'Rhabditina', 'Rhabditidae'),
    'Koerneria': ('Rhabditida', 'Diplogasteromorpha', 'Diplogastridae'),
    'Cooperia': ('Rhabditida', 'Strongylida', 'Cooperiidae'),
    'Dictyocaulus': ('Rhabditida', 'Strongylida', 'Dictyocaulidae'),
    'Heligmosomum': ('Rhabditida', 'Strongylida', 'Heligmosomidae'),
    'Nematodirus': ('Rhabditida', 'Strongylida', 'Molineidae'),
    'Oesophagostomum': ('Rhabditida', 'Strongylida', 'Chabertiidae'),
    'Ostertagia': ('Rhabditida', 'Strongylida', 'Haemonchidae'),
    'Trichostrongylus': ('Rhabditida', 'Strongylida', 'Trichostrongylidae'),
    'Acanthocheilonema': ('Rhabditida', 'Spirurina', 'Onchocercidae'),
    'Dirofilaria': ('Rhabditida', 'Spirurina', 'Onchocercidae'),
    'Litomosoides': ('Rhabditida', 'Spirurina', 'Onchocercidae'),
    'Spirocerca': ('Rhabditida', 'Spirurina', 'Spirocercidae'),
    'Thelazia': ('Rhabditida', 'Spirurina', 'Thelaziidae'),
    'Acrobeloides': ('Rhabditida', 'Tylenchina', 'Cephalobidae'),
    'Panagrellus': ('Rhabditida', 'Tylenchina', 'Panagrolaimidae'),
    'Panagrolaimus': ('Rhabditida', 'Tylenchina', 'Panagrolaimidae'),
    'Propanagrolaimus': ('Rhabditida', 'Tylenchina', 'Panagrolaimidae'),
    'Turbatrix': ('Rhabditida', 'Tylenchina', 'Panagrolaimidae'),
    'Sphaerularia': ('Rhabditida', 'Tylenchina', 'Sphaerulariidae'),
    'Daptonema': ('Monhysterida', 'Monhysterina', 'Xyalidae'),
    'Sphaerolaimus': ('Monhysterida', 'Monhysterina', 'Sphaerolaimidae'),
    'Echinomermella': ('Mermithida', 'Mermithina', 'Mermithidae'),
}

assemblies = []
for sp in species_dirs:
    sp_path = os.path.join(base, sp)
    analysis_dir = os.path.join(sp_path, 'analysis')
    if not os.path.exists(analysis_dir):
        continue
    try:
        entries = os.listdir(analysis_dir)
    except Exception:
        continue
    for acc in [e for e in entries if e.startswith('GCA_')]:
        acc_dir = os.path.join(analysis_dir, acc)
        busco_glob = glob.glob(os.path.join(acc_dir, 'pipelines/busco/*/results/busco'))
        if not busco_glob:
            continue
        busco_dir = busco_glob[0]
        chrom12_dir = os.path.join(busco_dir, 'chromadorea_odb12')
        if not os.path.exists(chrom12_dir):
            continue
        
        name_map_path = os.path.join(sp_path, 'assembly/release', acc, 'insdc', f'{acc}.name_mapping.tsv')
        num_chrs = 0
        chr_names = []
        if os.path.exists(name_map_path):
            with open(name_map_path) as f:
                for line in f:
                    p = line.strip().split('\t')
                    if len(p) >= 3 and p[2] and p[2] != 'MT':
                        num_chrs += 1
                        chr_names.append(p[1])
        
        if num_chrs == 0:
            continue
        
        # Get taxid from speciesops.json
        taxid = "NA"
        speciesops_path = os.path.join(sp_path, 'speciesops.json')
        if os.path.exists(speciesops_path):
            try:
                with open(speciesops_path) as f:
                    s_data = json.load(f)
                    taxid = str(s_data.get('taxon_id', 'NA'))
            except Exception:
                pass
        
        # Get BUSCO metrics from short_summary.json
        c_sc = 0
        c_du = 0
        c_fr = 0
        c_ms = 0
        summary_path = os.path.join(chrom12_dir, f'{acc}.chromadorea_odb12.short_summary.json')
        if os.path.exists(summary_path):
            try:
                with open(summary_path) as f:
                    b_data = json.load(f)
                    res = b_data.get('results', {})
                    c_sc = res.get('Single copy BUSCOs', 0)
                    c_du = res.get('Multi copy BUSCOs', 0)
                    c_fr = res.get('Fragmented BUSCOs', 0)
                    c_ms = res.get('Missing BUSCOs', 0)
            except Exception:
                pass
        
        genus = sp.split('_')[0]
        tax_info = taxonomy_lookup.get(genus, ('Nematoda', 'Unknown', 'Unknown'))
        
        assemblies.append({
            'species': sp,
            'accession': acc,
            'taxid': taxid,
            'order': tax_info[0],
            'suborder': tax_info[1],
            'family': tax_info[2],
            'n_chromosomes': num_chrs,
            'chr_names': chr_names,
            'busco_complete_sc': c_sc,
            'busco_duplicated': c_du,
            'busco_fragmented': c_fr,
            'busco_missing': c_ms,
            'busco_complete_pct': round((c_sc + c_du) / 1126.0 * 100, 1),
            'busco_sc_pct': round(c_sc / 1126.0 * 100, 1),
            'full_table_path': os.path.join(chrom12_dir, f'{acc}.chromadorea_odb12.full_table.tsv.gz'),
            'name_mapping_path': name_map_path
        })

print(f"Total Assemblies Processed: {len(assemblies)}")
# Output TSV header
print("Species\tAccession\tTaxID\tOrder\tSuborder\tFamily\tn_chrs\tChr_Names\tBUSCO_Complete_SC\tBUSCO_SC_Pct\tBUSCO_Total_Comp_Pct\tFull_Table_Path")
for a in sorted(assemblies, key=lambda x: (x['order'], x['suborder'], x['species'])):
    chrs_str = ";".join(a['chr_names'])
    print(f"{a['species']}\t{a['accession']}\t{a['taxid']}\t{a['order']}\t{a['suborder']}\t{a['family']}\t{a['n_chromosomes']}\t{chrs_str}\t{a['busco_complete_sc']}\t{a['busco_sc_pct']}%\t{a['busco_complete_pct']}%\t{a['full_table_path']}")
