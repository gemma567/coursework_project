Identify the most similar sequence

This Python script takes a mystery DNA sequence and aligns it with a database
of dog breed DNA sequences to find the highest percentage match. It also
generates a Pandas DataFrame with the percentage match for each dog breed
in the database and gives the p-value.


Requirements
This script requires the following packages to be installed:

Biopython
pandas
scipy
flask


Usage
To use this script, run breed_finder.py, providing your name and the path of the
file containing the mystery sequence to be tested as arguments.

Alternatively - for a better visual output - you can run app.py to see the results of as a webpage.
This serves the results of the mystery.fa alignment with the database at localhost:5000


Output
The script will generate a results folder with a file named your_name_results.txt.
This file will contain the following information:

The best match dog breed with its percentage match to the mystery sequence.
The p-value of the independent samples t-test.
A Pandas DataFrame with the percentage match for each dog breed in the database.

The Flask app will return a webpage containing:

The best match dog breed with its percentage match to the mystery sequence.
The p-value of the independent samples t-test.
A table showing the percentage match for each dog breed in the database.


Contributors
Gemma White


Contact Information
gwhite11@student.bbk.ac.uk