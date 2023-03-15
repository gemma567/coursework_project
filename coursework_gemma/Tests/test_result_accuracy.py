from Bio import SeqIO
from Bio.Align import PairwiseAligner


# creates a test sequence from a fasta file.
def read_genome(filename):
    test_genome = ""
    with open(filename, "r") as f:
        for line in f:
            if not line[0] == ">":
                test_genome += line.rstrip()
    return test_genome


# finds the best match and outputs the sequence with the best percentage
# match as well as the test sequence. These are then used to test
# the result with a different matching algorithm, using a database
# with the best match removed.
def best_match_test(test_seq, database_seq):
    """ a simple function to read in the test sequence and align
    it with the sequences in the database, finding the best match.
    It outputs the best matching sequence and the test sequence.
    These sequences are then passed to the align function to
    double-check the edit distance.
    """
    # Read in the mystery sequence of dog DNA
    test_seq = read_genome(test_seq)

    # Create a PairwiseAligner object
    aligner = PairwiseAligner()

    # Set the parameters for the aligner
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    # Create variables to keep track of the best match
    best_score = 0
    best_seq = ""
    # Loop through the dog breeds file and find the best match
    for record in SeqIO.parse(database_seq, 'fasta'):
        alignment = aligner.align(record.seq, test_seq)[0]  # get the first (best) alignment
        identity = alignment.score
        if identity > best_score:
            best_score = identity
            best_seq = record.seq

    # returns the closest matching sequence
    return best_seq, test_seq


mystery_file = "../Data_files/mystery.fa"
database = "../Data_files/dog_breeds.fa"
test_database = "../Data_files/dog_breeds_test_result.fa"

spaniel = best_match_test(mystery_file, database)
non_spaniel = best_match_test(mystery_file, test_database)


def align(x, y):
    """ A simple function that uses the Wagner Fisher algorithm to
    find the Levenshtein distance between the two strings. x and y
    will be the test sequence and the closest matching sequence
    returned from the best match test function.
    """
    if len(x) < len(y):
        x, y = y, x

    prev = list(range(len(y) + 1))
    curr = [0] * (len(y) + 1)

    for i, xc in enumerate(x, 1):
        curr[0] = i
        for j, yc in enumerate(y, 1):
            if xc == yc:
                curr[j] = prev[j-1]
            else:
                curr[j] = min(prev[j], curr[j-1], prev[j-1]) + 1
        prev, curr = curr, prev

    return prev[-1]


diffs = align(spaniel[0], spaniel[1])

print(f"There is an edit distance of {diffs} between the given sequence and the best match.")

diffs_2 = align(non_spaniel[0], non_spaniel[1])

print(f"There is an edit distance of {diffs_2} between the given sequence and the second best match.")

if diffs < diffs_2:
    print("Test passed.")
else:
    print("Test failed.")
