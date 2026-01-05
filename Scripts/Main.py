from pathlib import Path
import random

# ---------------- paths(best way I got this to work on both my ipad and PC, it was relative before) ----------------
base = Path(__file__).resolve().parent
dataDir = (base.parent / 'Data').resolve()
outDir = (base.parent / 'Results').resolve()
outDir.mkdir(exist_ok=True)

# ---------------- read input files ----------------
pro = open(dataDir / 'Proteins.fa', 'r')
codon_AA = open(dataDir / 'codon_AA.txt', 'r')

ProSeqs = []          # protein sequences
codonDict = {}        # amino acid -> list of codons

# ---------- read protein file ----------
seq = ""
for line in pro:
    if line.startswith(">"):
        if seq:
            ProSeqs.append(seq)
            seq = ""
        continue
    seq += line.strip()

if seq != "" and seq not in ProSeqs:                       # append last sequence
    ProSeqs.append(seq)
print(len(ProSeqs), "protein sequences loaded.")
choose = int(input(f"Choose protein index (0-{len(ProSeqs)-1}): "))      #im juust gonna trust the user here to put a valid index like cmon
pro.close()

# ---------- build codon dictionary ----------
for line in codon_AA:
    parts = line.strip().split()
    if len(parts) == 2:
        codon, amino_acid = parts
        if amino_acid not in codonDict:
            codonDict[amino_acid] = []
        codonDict[amino_acid].append(codon)

codon_AA.close()

# ---------------- funcs ----------------
#sequence from protein to rna (not DNA cause the file I had in trinket coddon_AA was RNA codons so wtv)
def backTranslate(proSeq):
    dnaSeq = ""
    for aminoAcid in proSeq:
        if aminoAcid in codonDict:
            dnaSeq += random.choice(codonDict[aminoAcid])
        else:
            dnaSeq += "NNN"
    return dnaSeq

#for the under 4 in a row condition thing
def hasLongRepeats(dna, maxRepeat=4):
    count = 1
    for i in range(1, len(dna)):
        if dna[i] == dna[i - 1]:
            count += 1
            if count > maxRepeat:
                return True
        else:
            count = 1
    return False

def gcContent(dna):
    if len(dna) == 0:   #just cause I had a division by 0 error once (I fixed the real issue but just in case)
        return 0
    return 100 * (dna.count('G') + dna.count('C')) / len(dna)

#between 40 and 60% gc cont
def isValidDNA(dna):
    gc = gcContent(dna)
    if gc < 40 or gc > 60:
        return False
    if hasLongRepeats(dna):
        return False
    return True

def RNA2DNA(rna):
  map = {'A': 'A', 'U': 'T', 'G': 'G', 'C': 'C'}
  dna = rna.upper()
  return ''.join(map.get(n, '') for n in dna)

#currently generates 5, can make n anything
def generateValidDNAs(proSeq, n=5):
    validSeqs = []
    attempts = 0

    while len(validSeqs) < n and attempts < 50000: #limit to not make it try to find valid forever
        dna = backTranslate(proSeq)
        if isValidDNA(dna):
            validSeqs.append((dna, gcContent(dna)))
        attempts += 1
    
    for i in range(len(validSeqs)):     #conver to DNA, no real reason to do this here but I just saw this and was like yeah thats good so peak
        rna_seq = validSeqs[i][0]
        dna_seq = RNA2DNA(rna_seq)
        validSeqs[i] = (dna_seq, validSeqs[i][1])

    return validSeqs

#finalized output thing
def outPutFasta(path, protein, dnaSeqs):
    with open(path, 'w') as f:
        f.write(">protein\n")
        f.write(protein + "\n")

        for i, (dna, gc) in enumerate(dnaSeqs, start=1):
            f.write(f">dna_{i} | GC={gc:.2f}%\n")
            f.write(dna + "\n")

# ---------------- main ----------------
protein = ProSeqs[choose]               #user chosen protein cause file has a bunch
dnaSeqs = generateValidDNAs(protein, n=5)

if len(dnaSeqs) < 5:
    print("warning: could not generate 5 valid sequences")      
else:
    outputFile = outDir / f"output(proSeq{choose+1}).fa"
    outPutFasta(outputFile, protein, dnaSeqs)
    print("output written to:", outputFile)
