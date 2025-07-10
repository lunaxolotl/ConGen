def filter_vcf_simple(input_file, output_file):
    samples_to_remove = {
    "valsusa_01", "rocciamelone_1", "ticino_1", "graubunden_1", "graubunden_2",
    "fribourg_1", "a.marittime_17", "a.marittime_18", "a.marittime_19",
    "a.marittime_20", "a.marittime_21", "lanzo_38", "lanzo_39",
    "granparadiso_46", "granparadiso_47", "granparadiso_48",
    "granparadiso_49", "granparadiso_50",
    # Added samples to remove:
    "AM_H", "APAMgoat1", "APAMgoat2", "APAMgoat3", "APAMgoat4",
    "Balme1t", "Balme2b",
    "VDG_1", "VDG_2", "VDG_3", "VDG_4", "VDG_5",
    "GPHB1", "HB_Lanzo_2021", "TI_ib",
    "GR_ib1", "GR_ib2", "FR_blanc"
}


    markers_to_remove = {
        (12, 20315678), (12, 40316282), (12, 70322246),
        (13, 10048151), (13, 25050600), (13, 55051511), (13, 80055158),
        (14, 25021789), (15, 20007044), (15, 65019388),
        (16, 30024980), (16, 50040434), (17, 30011647), (17, 45031043),
        (18, 15003075), (18, 55007685), (19, 10023896), (19, 45029681),
        (20, 30009781), (21, 20010967), (21, 55013563),
        (23, 25083623), (24, 20020388), (26, 40126121),
        (27, 20015442), (28, 5008953), (28, 40015753),
        (29, 20004050), (29, 35008446)
    }

    with open(input_file, 'rt') as fin, open(output_file, 'wt') as fout:
        total = kept = 0
        sample_indexes_to_keep = None

        for line in fin:
            if line.startswith('##'):
                # Metadata line, just copy
                fout.write(line)
                continue

            if line.startswith('#CHROM'):
                # Header line with sample names
                fields = line.strip().split('\t')
                fixed_cols = fields[:9]
                samples = fields[9:]

                # Determine indexes of samples to keep
                sample_indexes_to_keep = [i for i, s in enumerate(samples) if s not in samples_to_remove]

                # Write updated header line with filtered samples
                filtered_samples = [samples[i] for i in sample_indexes_to_keep]
                fout.write('\t'.join(fixed_cols + filtered_samples) + '\n')
                continue

            # Data line
            total += 1
            fields = line.strip().split('\t')
            chrom_raw = fields[0].replace('chr', '')
            pos = int(fields[1])

            # Convert chrom to int if possible
            try:
                chrom = int(chrom_raw)
            except ValueError:
                chrom = chrom_raw

            if (chrom, pos) in markers_to_remove:
                continue  # Skip this variant

            # Filter samples' genotype columns
            fixed_cols = fields[:9]
            sample_fields = fields[9:]
            filtered_samples = [sample_fields[i] for i in sample_indexes_to_keep]

            fout.write('\t'.join(fixed_cols + filtered_samples) + '\n')
            kept += 1

            if total % 10000 == 0:
                print(f"Processed {total:,} records, kept {kept:,}...")

    print(f"\nFiltering complete!")
    print(f"Total records processed: {total:,}")
    print(f"Records kept: {kept:,} ({kept/total:.1%})")
    print(f"Records removed: {total - kept:,} ({(total-kept)/total:.1%})")
    print(f"Filtered VCF saved to: {output_file}")


# Example usage:
input_vcf = "targets.dp20gq50.recode.vcf"  # or .vcf.gz (see note below)
output_vcf = "targets.dp20gq50.filtered.vcf"

filter_vcf_simple(input_vcf, output_vcf)
