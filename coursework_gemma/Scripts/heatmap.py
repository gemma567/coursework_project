from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
import seaborn as sns
import matplotlib.pyplot as plt

align = AlignIO.read("../Data_files/dog_msa.fas", "fasta")

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(align)

constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)

# create a heatmap of the distance matrix
sns.heatmap(dm, cmap='YlGnBu')
plt.title("Dog MSA Distance Matrix Heatmap")
plt.xlabel("Sequences")
plt.ylabel("Sequences")

# save the plot
plt.savefig("../results/heatmap.png", dpi=300)
