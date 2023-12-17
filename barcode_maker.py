from seqwalk import design
import pandas as pd
from Levenshtein import hamming
import sys
import pyfastx
from Bio.SeqUtils import MeltingTemp as mt
# python barcode_maker.py NUM_SEQ DISTANCE ITERATIONS > out.txt
#blacklist = set(['CG', 'GT', 'GG', 'GC'])
#pref = set(["AT", "TA", "GA", "AG"])

h = 'Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv'
m = 'Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv'
human = pd.read_csv(h, comment="#")
mouse = pd.read_csv(m, comment="#")
human = human['probe_seq'].str[5:-5].tolist()
mouse = mouse['probe_seq'].str[5:-5].tolist()

hv = 'Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv'
mv = 'Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv'
human_visium = pd.read_csv(hv, comment="#")
mouse_visium = pd.read_csv(mv, comment="#")
human_visium = human_visium['probe_seq'].str[5:-5].tolist()
mouse_visium = mouse_visium['probe_seq'].str[5:-5].tolist()

xenium = []
for name, seq in pyfastx.Fasta('xenium_human_brain_gene_expression_panel_v1_probe_sequences.fasta', build_index=False, uppercase=True):
    xenium.append(seq)
for name, seq in pyfastx.Fasta('xenium_mouse_brain_gene_expression_panel_v1.1_probe_sequences.fasta', build_index=False, uppercase=True):
    xenium.append(seq)
for name, seq in pyfastx.Fasta('xenium_mouse_multi_gene_expression_panel_v1_probe_sequences.fasta', build_index=False, uppercase=True):
    xenium.append(seq)
for name, seq in pyfastx.Fasta('xenium_human_multi_gene_expression_panel_v1_probe_sequences.fasta', build_index=False, uppercase=True):
    xenium.append(seq)

merged = human + mouse + human_visium + mouse_visium + xenium
merged = set(merged)
left = [w[:20] for w in merged]
right = [w[:-20] for w in merged]

with open("re.txt") as file:
    res = [line.rstrip() for line in file]

library = set()
for it in range(int(sys.argv[3])):
    library.update(design.max_orthogonality(int(sys.argv[1]), 40, alphabet="ACGT", RCfree=True, GClims=(15, 35), prevented_patterns=res))

dist = int(sys.argv[2])
j = 0
for i in library:
    melt = mt.Tm_NN(i)
    lig = i[19:21]
    if lig != "AT" or lig != "TA" or lig != "GA" or lig != "AG" or melt < 68 or melt > 82:
        continue
    flag = False
    for probe in left:
        l = i[:20]
        ham = hamming(l, probe)
        melt = mt.Tm_NN(l)
        if ham < dist or melt < 50 or melt > 70:
            flag = True
            break
    if flag:        
        continue
    for probe in right:
        r = i[:-20]
        ham = hamming(r, probe)
        melt = mt.Tm_NN(r)
        if ham < dist or melt < 50 or melt > 70:
            flag = True
            break
    if flag:
        continue
    print(i)
    j += 1
    if j == 201:
        break
