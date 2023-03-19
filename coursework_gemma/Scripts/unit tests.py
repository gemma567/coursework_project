import unittest
import os
from dog_module import align_dog_breeds, read_genome


class TestAlignDogBreeds(unittest.TestCase):
    def test_align_dog_breeds(self):
        test_sequence = "../Data_files/test_sequence.txt"
        database_seqs = "../Data_files/test_database.fa"
        result = align_dog_breeds(test_sequence, database_seqs)
        self.assertEqual(result[0], "Portuguese Pointing Dog-long hair")
        self.assertAlmostEqual(result[1], 100.00, places=2)
        self.assertLess(result[2], 1)
        self.assertIsInstance(result[3], dict)


class TestReadGenome(unittest.TestCase):

    def test_read_genome(self):
        # Defines a test input file with known content.
        test_file = 'test_file.txt'
        with open(test_file, 'w') as f:
            f.write('>Header1\n')
            f.write('ATCG\n')
            f.write('>Header2\n')
            f.write('CGTA\n')

        # Calls the function using the test file.
        result = read_genome(test_file)

        # Checks the output is as expected.
        expected_output = 'ATCGCGTA'
        self.assertEqual(result, expected_output)

        # Deletes the test file.
        os.remove(test_file)


if __name__ == '__main__':
    unittest.main()

