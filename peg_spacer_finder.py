import numpy as np
import pandas as pd
import re

# Using a raw string to ensure backslashes are interpreted correctly
data = pd.read_csv(r"PrimeDesign_Pooled_20241031_07.46.21.51.csv")

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
spacer_sequences = data["Spacer_sequence"]
pams = data["PAM_sequence"]
flipped_sequences = [] # tuple of flipped reference and alternate sequences
# need to flip some of the strands

# NEED TO make cleaned_sequences contain reference_sequence, alternate_sequence
# let's just make a reference and alternate for every strand

#reference is keep left of / and delete the right


reference_strands = []
for sequence in original_sequences: # THIS WORK
    temp = re.sub(r'/(.*?)\)', '', sequence) # hopefully removing just the right
    temp = temp.replace( '(', '') # removing )
    reference_strands.append(temp)

alternate_strands = []
for sequence in original_sequences: # THIS WORKS
    temp = re.sub(r'\((.*?)/', '', sequence) # hopefully removing just the left
    temp = temp.replace( ')', '') # removing )
    alternate_strands.append(temp)

# now we have all the reference and alternate strands, now we need the spam
spacer_pams = [f"{s}{p}" for s, p in zip(spacer_sequences, pams)]  # Combined spacer and PAM WORKS

# we have reference sequences, alternate sequence, spams, and peg extensions

flipper = str.maketrans("ACTG", "TGAC","()")
for reference, alternate in zip(reference_strands, alternate_strands): # now we have the flipped strands
    t1 = reference.translate(flipper)[::-1]
    t2 = alternate.translate(flipper)[::-1]
    flipped_sequences.append((t1, t2))





pattern = r'/.*?\)'
cleaned_sequences = [re.sub(pattern, '', seq) for seq in original_sequences]

'''for s,p in zip(spacer_sequences,pams):
    mixed = spacer_sequences + pams
    spacer_pams.append(mixed) # spacer_pams has spacers + pams'''

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
#print(locations)





# spam matches to the reference (keep the left and delete the right)
# peg extstension matches to the alternate sequence (keep the right and delete the left)
    # with peg extensions, since it matches the alternate, the final coordinate you just add on the difference between the original and the edit (og length - edit length)


''' need to create an array of both the reference and the alternate sequences, that way I can just easily compare 
    spam --> reference and peg extension --> alternate'''