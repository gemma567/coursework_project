from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo import draw
import matplotlib.pyplot as plt

align = AlignIO.read("../Data_files/dog_msa.fas", "fasta")

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(align)

constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)


# define a function to label the nodes
def label_func(clade):
    if clade.name:
        return clade.name
    else:
        return ""


plt.figure(figsize=(10, 10))

# draw the tree and customize it
draw(tree, label_func=label_func, do_show=False, show_confidence=False)
plt.title("Dog MSA Tree")
plt.xlabel("Distance")
plt.ylabel("Sequences")

# save the plot
plt.savefig("../results/tree.png", dpi=300)


