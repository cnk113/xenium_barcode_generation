x = readLines("library.txt")

idx = paste0("AATCTAGA",x[1:2],"GACGTCAA")

bit1 = paste0("AAACGCGT",x[3:42],"TCCGGAAA")

bit2 = paste0("AAGGATCC",x[43:82],"GTCGACAA")

bit3 = paste0("AAGGTACC",x[83:122],"TGATCAAA")

bit4 = paste0("AACGTACG",x[123:162],"CCCGGGAA")

bit5 = paste0("AAATGCAT",x[163:202],"CCTGCAGGAA")

z = c(idx, bit1, bit2, bit3, bit4, bit5)

writeLines(z, "library_idt.txt")
