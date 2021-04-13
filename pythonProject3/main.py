from DNAtoolkit import *
from utilities import colored
import random
rndDNA = "".join([random.choice(Nucleotides)
                  for nuc in range(1000)])

DNAStr = validateSeq(rndDNA)
print(f'\nSequence: {colored(DNAStr)}\n')

print(f'[1] Sequence Length: {len(DNAStr)}\n')

print(colored(f'[2] Nucleotide Frequency: {countNucFrequency(DNAStr)}\n'))

print(f'[3] Transcription: {colored(transcription(DNAStr))}\n')

# 3 lines of print together
print(f"[4] DNA String + Reverse Complement:\n5' {colored(DNAStr)} 3' ")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {colored(reverse_complement(DNAStr)[::-1])} 5'\n")

print(f'[5] GC content: {gc_content(DNAStr)}%\n')
print(f'[6] GC content in Subsection k=5: {gc_content_subset(DNAStr, k=50)}%\n')
print(f'[7] Aminoacids Sequence from DNA: {translate_seq(DNAStr, 1)}\n')
print(f'[8] Codon Frequency: {codon_usage(DNAStr, "_")}\n')
print(f'[9] Reading frames:')
for frame in gen_reading_frames(DNAStr):
    print(frame)

print('\n[10] All prots in 6 open reading frames:')
for prot in all_prot_from_orf(DNAStr, 0, 0, True):
    print(f'{prot}')


