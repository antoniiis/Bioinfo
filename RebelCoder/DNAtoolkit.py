# DNA Toolkit file
from collections import *
from structures import *

# Check the sequence to make sure it is a DNA string
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
       tmpFreqDict[nuc] += 1
    return tmpFreqDict
    # Or return dict(collections.Counter(seq))


def transcription(seq):
    # turns DNA into RNA
    return seq.replace("T", "U")


def reverse_complement(seq):
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]
    #Second way
    #return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]


def gc_content(seq):
    #GC content in a DNA/RNA sequence
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)


def gc_content_subset(seq, k=20):
    #GC content of a DNA/RNA sub-sequence length k... k=20 by default
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i +k]
        res.append(gc_content(subseq))
    return res


def translate_seq(seq, init_pos= 0):
    #Translates a DNA sequence into an aminoacid sequence
    return [DNA_Codons[seq[pos:pos +3]] for pos in range(init_pos, len(seq) -2, 3)]


def codon_usage(seq, aminoacid):
    #Provides the frequency of each codon encoding a given aminoacid in a DNA sequence
    tmplist = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i +3]] == aminoacid:
            tmplist.append(seq[i:i +3])

    freqDict = dict(Counter(tmplist))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict


def gen_reading_frames(seq):
    #Generate the six reading frames of a DNA sequence, including the reverse complement
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames

def proteins_from_rf(aa_seq):
    #Compute all possible proteins in an aminoacid seq and return a list of possible proteins
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            #STOP accumulating aminoacids if _ was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            #START accumulating aminoacids if M was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_prot_from_orf(seq, startReadPos= 0, endReadPos= 0, ordered= False):
    #Compute all possible proteins for all open reading frames
    #Protein search DB:
    #API can be used to pull protien info
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = gen_reading_frames(seq)

    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)

    if ordered:
        return sorted(res, key= len, reverse= True)
    return res


