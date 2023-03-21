from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo import draw
import matplotlib.pyplot as plt

align = AlignIO.read("../Data_files/dog_msa.fas", "fasta")
# print(align)


calculator = DistanceCalculator('identity')
dm = calculator.get_distance(align)


constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)

plt.figure(figsize=(10, 10))

# Draw the tree and save to a file
draw(tree, do_show=False, show_confidence=False)
plt.savefig("../results/tree.png", dpi=300)
