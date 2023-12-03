from seqwalk import design
import pandas as pd
from Levenshtein import hamming
import sys
import pyfastx

# python barcode_maker.py NUM_SEQ DISTANCE ITERATIONS > out.txt
#blacklist = set(['CG', 'GT', 'GG', 'GC'])
#pref = set(["AT", "TA", "GA", "AG"])

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

library = set()
for it in range(int(sys.argv[3])):
    library.update(design.max_orthogonality(int(sys.argv[1]), 40, alphabet="ACGT", RCfree=True, GClims=(16, 29),
        prevented_patterns=["AAAA", "CCCC", "GGGG","TTTT","ACGCGT","ATGCAT","GGATCC","GTCGAC","GGTACC","TGATCA","CGTACG","GCGCGC","ATGCAT","CCTGCAGG","TCTAGA","GACGTC"]))
lib = []
dist = int(sys.argv[2])
for i in library:
    flag = False
    if i[19:21] == "AT" or i[19:21] == "TA" or i[19:21] == "GA" or i[19:21] == "AG":
        for probe in left:
            l = i[:20]
            ham = hamming(l, probe)
            gc = l.count("G") + l.count("C")
            if ham > dist and gc > 7 and gc < 16:
                flag = True
                break
        if flag:
            continue
        for probe in right:
            r = i[:-20]
            ham = hamming(r, probe)
            gc = r.count("G") + r.count("C")
            if ham > dist and gc > 7 and gc < 16:
                 flag = True
                 break
        if flag:
            continue
    print(i)
