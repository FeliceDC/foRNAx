#!/usr/bin/env python3
import sys
import pandas as pd
import immucellai2

fc_file = sys.argv[1]
threads = sys.argv[2] if len(sys.argv) > 2 else 4
# Cattura il parametro del tessuto da Nextflow (default: tumor)
tissue = sys.argv[3].lower() if len(sys.argv) > 3 else "tumor"

print("1. Caricamento della matrice featureCounts...")
df = pd.read_csv(fc_file, sep='\t', comment='#', index_col=0)

lengths = df['Length']
counts = df.iloc[:, 5:]

counts.columns = [c.replace('.Aligned.sortedByCoord.out.bam', '').replace('.bam', '').lstrip('X') for c in counts.columns]

print("2. Calcolo dei TPM (Transcripts Per Million)...")
rpk = counts.div(lengths / 1000, axis=0)
tpm = rpk.div(rpk.sum(axis=0) / 1e6, axis=1)

tpm_file = "tpm_matrix.txt"
tpm.to_csv(tpm_file, sep='\t')

print(f"3. Avvio di ImmuCellAI 2.0 (Deconvoluzione con reference: {tissue.upper()})...")

if tissue == "blood":
    ref_data = immucellai2.load_blood_reference_data()
else:
    ref_data = immucellai2.load_tumor_reference_data()

immucellai2.run_ImmuCellAI2(
    reference_file=ref_data,
    sample_file=tpm_file,
    output_file="ImmuCellAI2_results.xlsx",
    thread_num=int(threads)
)
print("Analisi completata con successo!")
