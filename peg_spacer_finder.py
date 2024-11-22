import numpy as np
import pandas as pd
import re

# Using a raw string to ensure backslashes are interpreted correctly
data = pd.read_csv(r"C:\Users\that9\OneDrive\Documents\PrimeDesign_Pooled_20241031_07.46.21.51.csv")

# data has all that data 
# for now let's just try to find peg extensions in the original sequence

# og sequence column = target sequence COLUMN 3 / C
# peg extension column = pegRNA_extension COLUMN 12 / L
# make array of all target sequence column

original_sequences = data["Target_sequence"]
flips = data["Strand"]
adjusted_sequences = []
peg_extensions = data["pegRNA_extension"]
locations = [] # add tuples of locations
# need to flip some of the strands

pattern = r'/.*?\)'
cleaned_sequences = [re.sub(pattern, '', seq) for seq in original_sequences]
spacer_sequences = data["Spacer_sequence"]
pams = data["PAM_sequence"]
# spacer_pams = []
spacer_pams = [f"{s}{p}" for s, p in zip(spacer_sequences, pams)]  # Combined spacer and PAM
'''for s,p in zip(spacer_sequences,pams):
    mixed = spacer_sequences + pams
    spacer_pams.append(mixed) # spacer_pams has spacers + pams'''
flipper = str.maketrans("ACTG", "TGAC","()")
for sequence, sign in zip(cleaned_sequences, flips): # pairs each sequence with its sign
    if (sign == "-"): # need to make the strand antiparallel and complementary 
        adjusted_sequences.append(sequence.translate(flipper)[::-1])
    else: # +, so use normal sequence
        adjusted_sequences.append(sequence) # use this now

# adjusted sequences now has all the necessary flips
# + is reverse for extensions, - is normal for extension

for flip, clean, sequence, peg, spacer in zip(flips, cleaned_sequences, adjusted_sequences, peg_extensions, spacer_pams): # need to find where the substring is
    if (flip == '+'): #use reverse compliment
        index = sequence.find(peg)
    else: index = clean.find(peg)
    index2 = sequence.find(spacer)
    #print(index)

    lo = []

    if (index != - 1): # BINGO -- have the initial index of the spot
        lo.append(index)
        lo.append(index + len(sequence))
    else: # not in the sequence
        lo.append("X")
        lo.append("X")
    if (index2 != -1): # BINGO but for spacerPAM
        lo.append(index2)
        lo.append(index2 + len(spacer))
    else:
        lo.append("X")
        lo.append("X")
    locations.append(lo)


locations_df = pd.DataFrame(locations, columns = ['Peg_Start', 'Peg_End', 'SPAM_Start', 'SPAM_End'])
output_file_path = r"C:\Users\that9\OneDrive\Documents\output_locations.csv"
locations_df.to_csv(output_file_path, index=False)
print(locations)
# bruh burh
