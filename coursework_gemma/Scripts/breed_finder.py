import os
from dog_module import align_dog_breeds
from dog_module import return_results

# gets the users name - this will be the name of the result file.
name = input("Please enter your name: ")

# gets filename from user - (please copy the full path of mystery.fa in the folder: Data_files to test my code!)
test_seq = input("Please enter full path of sequence filename: ")

# converts filename to Windows format if required.
if os.name == "nt":
    test_seq = test_seq.replace("/", "\\")

results = align_dog_breeds(test_seq)

my_alignment = return_results(results, name)
