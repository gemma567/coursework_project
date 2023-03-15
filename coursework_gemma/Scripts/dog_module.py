from Bio import SeqIO
from Bio.Align import PairwiseAligner
import random
import scipy.stats as stats
import os
import pandas as pd


def read_genome(filename):
    """ a function to read in a file containing a DNA sequence,
    strip it of any lines starting > and return just the sequence.

    input:
    filename: file containing sequence to be tested.

    output:
    a sequence stripped of lines containing >

    """
    test_genome = ""
    with open(filename, "r") as file:
        for line in file:
            if not line[0] == ">":
                test_genome += line.rstrip()
    return test_genome


def align_dog_breeds(test_sequence, database_seqs="../Data_files/dog_breeds.fa"):
    """ a function that reads in a mystery sequence and a database
    of sequences, aligns the sequences and finds the best match
    as well as printing out the percentage match of the mystery
    sequence with each sequence in the database and the p-value.

    input:
    test_seq: str - filename of sequence to be tested
    database_seqs: str - filename of database sequence.
        if this is left blank, the function will use the dog_breeds.fa
        file as a default.

    output:
    best_breed: str - breed of best match in database
    best_percentage_chance: float - percentage match of best match
    p_value: float - p value of best match
    doggy_dict: dict - a dictionary containing the percentage matches across
    the database.

    """
    test_sequence = read_genome(test_sequence)
    # Create a PairwiseAligner object
    aligner = PairwiseAligner()

    # Set the parameters for the aligner
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    # initialises variables and creates empty lists for the sequence percentage
    # identity and random sequence percentage identity for the t-test.
    percentages_list = []
    random_percentage_list = []
    doggy_dict = {}
    best_score = 0
    best_breed = ""
    best_percentage_chance = ""
    # iterates through the sequences in the file, finds the alignment score
    # and turns it into a percentage.
    # for percentage - len of sequence is * 2 since each match is plus 2
    # in the aligner parameters.
    for record in SeqIO.parse(database_seqs, 'fasta'):
        breed = record.description.split("] [")[6][6:]
        alignment = aligner.align(record.seq, test_sequence)[0]  # get the first (best) alignment
        identity = alignment.score
        if identity > best_score:
            best_score = identity
            best_breed = breed
            best_percentage_chance = (alignment.score / (len(record.seq) * 2)) * 100
        percentage_chance = round((alignment.score / (len(record.seq) * 2)) * 100, 2)
        percentages_list.append(percentage_chance)
        dog_id = record.id

        # generates a random sequence of the same length as the record
        random_seq = ''.join(random.choices(['A', 'C', 'G', 'T'], k=len(record.seq)))

        # calculates the percentage chance of a match with the random sequence
        alignment = aligner.align(random_seq, test_sequence)[0]  # get the first (best) alignment
        random_percentage_chance = (alignment.score / (len(random_seq) * 2)) * 100
        random_percentage_list.append(random_percentage_chance)

        # creates dictionary with the percentages across database
        doggy_dict[dog_id] = (breed, float(percentage_chance))
    # performs the independent samples t-test on the two sets of data we created in
    # the lists and that t-test is used to generate the p-value
    t_stat, p_value = stats.ttest_ind(random_percentage_list, percentages_list)
    # returns the results as a tuple
    return best_breed, best_percentage_chance, p_value, doggy_dict


def return_results(dog_align_result, name):
    """ A simple function that saves the results of the
     align_dog_breeds function to a file.
    it creates the results folder if it does not already exist
    and then writes the results of the align_dog_breeds
    function to a file named your_name_results.

    input:

    name: str - name of the user - used to create a unique filename.

    dog_align_results: - results of the align_dog_breeds function.

    """
    # gets the parent directory of the current directory
    parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

    # creates a folder called 'results' inside the parent directory if it doesn't already exist
    results_folder = os.path.join(parent_dir, 'results')
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    filename = name + '_results.txt'
    # sets the path for the results file inside the 'results' folder
    path = os.path.join(results_folder, filename)

    # opens the file at the given path for writing
    with open(path, 'w') as f:
        # writes the results to the file 'results.txt'
        f.write(f"Your results are as follows:")
        f.write('\n')
        f.write(f"The best match is {dog_align_result[0]}, with an identity of"
                f" {dog_align_result[1]:.2f}% and a p-value of {dog_align_result[2]}")

        # creates a Pandas DataFrame from the doggy_dict dictionary
        # this improves the readability of the output in the file
        df = pd.DataFrame.from_dict(dog_align_result[3], orient='index',
                                    columns=['Dog sequence id and Breed', 'Percentage match'])

        # writes the DataFrame to the file
        f.write('\n\nPercentage across database:')
        f.write('\n')
        df.to_csv(f, sep='\t')
