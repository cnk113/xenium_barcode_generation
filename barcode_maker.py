from seqwalk import design
import pandas as pd
from Levenshtein import hamming
import sys
import pyfastx

# python barcode_maker.py NUM_SEQ DISTANCE ITERATIONS > out.txt

#blacklist = set(['CG', 'GT', 'GG', 'GC'])
pref = set(["AT", "TA", "GA", "AG"])

h = 'Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv'
m = 'Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv'
human = pd.read_csv(h, comment="#")
mouse = pd.read_csv(m, comment="#")
human = human['probe_seq'].str[5:-5].tolist()
mouse = mouse['probe_seq'].str[5:-5].tolist()

xenium = []
for name, seq in pyfastx.Fasta('xenium_human_brain_gene_expression_panel_v1_probe_sequences.fasta', build_index=False, uppercase=True):
    xenium.append(seq)
for name, seq in pyfastx.Fasta('xenium_mouse_brain_gene_expression_panel_v1.1_probe_sequences.fasta', build_index=False, uppercase=True):
    xenium.append(seq)

merged = human + mouse + xenium

left = [w[:20] for w in merged]
right = [w[:-20] for w in merged]

dist = int(sys.argv[2])
lib = []
for it in range(int(sys.argv[3])-1):
    library = design.max_orthogonality(int(sys.argv[1]), 40, alphabet="ACGT", RCfree=True, GClims=(16, 28),
        prevented_patterns=["AAAAA", "CCCCC", "GGGGG","TTTTT","ACGCGT","ATGCAT","GGATCC","GTCGAC","GGTACC","TGATCA","CGTACG","GCGCGC","ATGCAT","CCTGCAGG","TCTAGA","GACGTC"])
    flag = False
    for i in library:
        if i[19:21] not in pref:      
            for probe in left:
                distance = hamming(i[:20], probe)
                if distance < dist:
                    flag = True
                    break
            if flag:
                continue
            flag = False
            for probe in right:
                distance = hamming(i[:-20], probe)
                if distance < dist:
                    flag = True
                    break
            if flag:
                continue
        lib.append(i)
for i in lib:
    print(i)

