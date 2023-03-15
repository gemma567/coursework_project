from Bio import SeqIO
from Bio.Align import PairwiseAligner
import random
import scipy.stats as stats

# This script tests the accuracy of the p-value in the dog_module function. For the mystery
# file and the dog_breed file, the dog_module returns a p-value of zero which means the p-value
# is so close to 0 that the computer is interpreting as a 0. I wanted to check this is correct by creating
# a list with zero percentages and checking against the result percentages - which would certainly return
# a result of 0 (or so close to 0 that 0 is returned). I also tested the script using the exact same values
# in each data set for which we would expect a return p-value of 1.


def read_genome(filename):
    test_genome = ""
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                test_genome += line.rstrip()
    return test_genome


# Read in the mystery sequence of dog DNA
test_seq = read_genome('../Data_files/mystery.fa')

# Create a PairwiseAligner object
aligner = PairwiseAligner()

# Set the parameters for the aligner
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5

percentages_list = []
random_percentage_list = []
percentage_list_zero = []
data_for_analysis = {}
# Loop through the dog breeds file and saves the percentage chance of
# a match in the percentages_list.
for record in SeqIO.parse('../Data_files/dog_breeds.fa', 'fasta'):
    breed = record.description.split("] [")[6][:]
    alignment = aligner.align(record.seq, test_seq)[0]
    percentage_chance = (alignment.score / (len(record.seq) * 2)) * 100
    percentages_list.append(percentage_chance)
    # Generates a random sequence of the same length as the record
    random_seq = ''.join(random.choices(['A', 'C', 'G', 'T'], k=len(record.seq)))
    # Calculates the percentage chance of a match with the random sequence
    # as well as a list of percentages that are all 0.
    alignment = aligner.align(random_seq, test_seq)[0]
    random_percentage_chance = (alignment.score / (len(random_seq) * 2)) * 100
    random_percentage_list.append(random_percentage_chance)
    percentage_list_zero.append(0)


# creates an identical copy of the list.
random_percentage_list_copy = random_percentage_list.copy()

# the following will test if the t-test is working by comparing the list
# with an identical list. Result should be 1.

t_stat_1, p_value_1 = stats.ttest_ind(random_percentage_list, random_percentage_list_copy)

# Print the results
print("t-statistic:", t_stat_1)
print("p-value:", p_value_1)

# this test compares the random list with the list of zeros. The result
# should be zero.

t_stat_2, p_value_2 = stats.ttest_ind(random_percentage_list, percentage_list_zero)

# Print the results
print("t-statistic:", t_stat_2)
print("p-value:", p_value_2)

if p_value_1 == 1 and p_value_2 == 0.0:
    print("Test passed")
else:
    print("Test failed")