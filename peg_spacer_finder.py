import numpy as np
import pandas as pd
import re
import random

# Using a raw string to ensure backslashes are interpreted correctly
data = pd.read_csv(r"PrimeDesign_Pooled_20241031_07.46.21.51.csv")
random.seed(18) # THIS IS THE SEED
# data has all that data 
# for now let's just try to find peg extensions in the original sequence

# og sequence column = target sequence COLUMN 3 / C
# peg extension column = pegRNA_extension COLUMN 12 / L
# make array of all target sequence column


# CREATING TEST CASE ISH ---------         
tester = "AAA(TTT/CCCCC)GGG"
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
for sequence in original_sequences: # THIS WORK!!
    temp = re.sub(r'/(.*?)\)', '', sequence) # hopefully removing just the right
    temp = temp.replace( '(', '') # removing )
    reference_strands.append(temp)

alternate_strands = []
for sequence in original_sequences: # THIS WORKS!!
    temp = re.sub(r'\((.*?)/', '', sequence) # hopefully removing just the left
    temp = temp.replace( ')', '') # removing )
    alternate_strands.append(temp)

# now we have all the reference and alternate strands, now we need the spam
spacer_pams = [f"{s}{p}" for s, p in zip(spacer_sequences, pams)]  # Combined spacer and PAM WORKS!!!
#****** ^^^ potential issue her with the reverse comp on the -????????? *******

# we have reference sequences, alternate sequence, spams, and peg extensions

flipper = str.maketrans("ACTG", "TGAC","()")
for reference, alternate in zip(reference_strands, alternate_strands): # now we have the flipped strands
    t1 = reference.translate(flipper)[::-1]
    t2 = alternate.translate(flipper)[::-1]
    flipped_sequences.append((t1, t2))

###### if the peg starting postiion is >201 or if the end is >201, then we need to change the indexes 
    # need to add (original length - edit length) to the index
    # for the alternate sequence, need to just take the (reference insert - edit insert) and add to array

length_differences = []


for sequence in original_sequences:
    # get original bases from ( to /
    match1 = re.search(r'\((.*?)/', sequence)
    match2 = re.search(r'/(.*?)\)', sequence)

    if match1:
        original_bases = match1.group(0)
        length1 = len(original_bases)
    if match2:
        alternate_bases = match2.group(0)
        length2 = len(alternate_bases)
    difference = length1 - length2
    length_differences.append(difference) # cause ( / / )



# if first term or last term is above 201, add the difference
count = 0
for reference, alternate, spam, peg, flip, difference in zip(reference_strands, alternate_strands, spacer_pams, peg_extensions, flips, length_differences): # need to find where the substring is
    lo = []
    # the spam is reverse compliment on -, the peg is reverse compliment on +
    if (flip == '+'): 
        index1 = reference.find(spam)
        #index2 = flipped_sequences[count][1].find(peg)
        index2 = (len(flipped_sequences[count][1]) - flipped_sequences[count][1].find(peg) - len(peg))
        
    elif (flip == '-'):
        index1 = (len(flipped_sequences[count][0]) - flipped_sequences[count][0].find(spam) - len(spam))
        index2 = alternate.find(peg)
    else:
        raise ("YOUR DATA HAS FORMATTING ERRORS")
    

    if (index1 != -1): # BINGO -- we have a match for the reference + spam
        lo.append(index1)
        lo.append(index1 + len(spam))
    else: # NO MATCH -- no match for reference + spam
        lo.append("X")
        lo.append("X")

    if (index2 != -1): # BINGO -- we have a match for the alternate + peg
        back = index2 + len(peg)
        if (index2 > 201):
            index2 = index2 + difference
        if (back > 201):
            back = back + difference
        lo.append(index2)
        lo.append(back)
    else: # NO MATCH -- no match for alternate + peg
        lo.append("X")
        lo.append("X")
    smallest = min(lo)
    biggest = max(lo)
    lo.append(smallest)
    lo.append(biggest)
    extracted = reference[smallest: (biggest)]
    if (flip == '-'):
        extracted = extracted.translate(flipper)[::-1]
    lo.append(spam)
    lo.append(peg)
    lo.append(extracted)
    locations.append(lo) # add the loactions to the location list
    count = count + 1



#print(locations)



#SPAM IS REVERSE COMP ON -
# PEG IS REVERSE COMP ON +



# need to eventually get shortest start and furthest end

# spam matches to the reference (keep the left and delete the right)
# peg extstension matches to the alternate sequence (keep the right and delete the left)
    # with peg extensions, since it matches the alternate, the final coordinate you just add on the difference between the original and the edit (og length - edit length)


''' need to create an array of both the reference and the alternate sequences, that way I can just easily compare 
    spam --> reference and peg extension --> alternate'''

# LFG the pseudo targets are properly created 
# now need to make the final oligo

''' final oligo should be u6_stub + protospacer (no PAM) + FE_hairpin + extension + 
tevO-Q1_motif_and_terminator + Reporter (pseudo target) + 6bp barcode + R2seq_stub'''

u6_stub = "TGTGGAAAGGACGAAACACCG"
FE_hairpin = "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
tev = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"
seq_stub = "AGATCGGAAGAGCACACGTCT"
oligos = []
used_codes = []
barMaker = str.maketrans("1234", "ACGT")
index = 0
# pseudo @ locations[count][8]
for count in locations:
    pseudo = count[8]
    # need to generate barcode
    while True:
        random_number = str(int(''.join(str(random.randint(1,4)) for num in range(6))))
        if (random_number not in used_codes): #FRESH CODE
            used_codes.append(random_number)  # add to used list
            break
    barcode = random_number.translate(barMaker)
    print(barcode)

    oligos.append(u6_stub + spacer_sequences[index] + FE_hairpin + peg_extensions[index] + tev + pseudo + barcode + seq_stub)
    index = index + 1




locations_df = pd.DataFrame(locations, columns = ['SPAM_Start', 'SPAM_End', 'Peg_Start', 'Peg_End', 'Earliest Start', 'Latest End',"S+PAM", "Peg Extension", "Pseudo Target"])
locations_df['Oligos'] = oligos
output_file_path = "output_locations.csv"
locations_df.to_csv(output_file_path, index=False)



