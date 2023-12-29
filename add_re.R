library(Biostrings)

x = readLines("library.txt")

idx = paste0("CTAGA",x[1:2],"GACGT")
dna = DNAStringSet(x[1:2])
idx_rc = paste0("C",reverseComplement(dna),"T")

bit1 = paste0("CGCGT",x[3:42],"T")
dna = DNAStringSet(x[3:42])
bit1_rc = paste0("CCGGA",reverseComplement(dna),"A")

bit2 = paste0("GATCC",x[43:82],"G")
dna = DNAStringSet(x[43:82])
bit2_rc = paste0("TCGAC",reverseComplement(dna),"G")

bit3 = paste0("C",x[83:122],"T")
dna = DNAStringSet(x[83:122])
bit3_rc = paste0("GATCA",reverseComplement(dna),"GGTAC")

bit4 = paste0("GTACG",x[123:162],"C")
dna = DNAStringSet(x[123:162])
bit4_rc = paste0("CCGGG",reverseComplement(dna),"C")

bit5 = paste0("T",x[163:202],"CCTGCA")
dna = DNAStringSet(x[163:202])
bit5_rc = paste0("GG",reverseComplement(dna),"ATGCA")

bc = list(bit1, bit2, bit3, bit4, bit5)
bc_rc = list(bit1_rc, bit2_rc, bit3_rc, bit4_rc, bit5_rc)
b = list("Bit1_", "Bit2_", "Bit3_", "Bit4_", "Bit5_")
plates = lapply(c(1:5), function(i){
	plate = expand.grid(c("A","B","C","D","E","F","G","H"), c(1:5))
        plate = rbind(plate, plate)
	plate$Well = paste0(plate$Var1, plate$Var2)
	name = expand.grid(c(1:40), c("A","B"))
	plate$Name = paste0(b[[i]], name$Var1, name$Var2)
	plate$Sequence = c(bc[[i]],bc_rc[[i]]) 
        write.csv(plate[,c(3:5)], paste0(b[[i]], "plate.csv"), row.names=FALSE)
        plate
})
idx_plate = data.frame(Well = c("A6","B6","A6","B6"),Name = c("Idx1A","Idx2A","Idx1B","Idx2B"))
idx_plate$Sequence = c(idx, idx_rc)
plate1 = rbind(plates[[1]][,c(3:5)], idx_plate)
write.csv(plate1, paste0(b[[1]], "plate.csv"), row.names=FALSE)


