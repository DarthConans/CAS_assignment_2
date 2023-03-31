import matplotlib.pyplot as plt
from matplotlib import style
from utils import *

sequence = load_sequence(False)
synonymous = []
nonsynonymous = []

xs = list(range(1, 11))
for hop in xs:
    syn, nonsyn = get_antigenic_synonymous_nonsynonymous_mutations(sequence, hop)
    synonymous.append(syn)
    nonsynonymous.append(nonsyn)

style.use('fivethirtyeight')

plt.rc('axes', titlesize=14)
plt.rc('axes', labelsize=12)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)

fig = plt.figure(figsize=(8, 6))
plot = fig.add_subplot(1, 1, 1)

plot.set_title("Antigenically Neutral vs Non-Antigenically Neutral Mutations\nin Neighborhood at n-hop")

plot.set_xlabel("n (steps)")
plot.set_xticks(range(0, int(max(xs))+1, 1))

plot.set_ylabel("unique mutations in neighborhood")
plot.set_yscale("log")

plt.tight_layout()

plot.plot(xs, synonymous, linewidth=1, marker="o", color="red", label="Antigenically Neutral Mutations")

plot.plot(xs, nonsynonymous, linewidth=1, marker="o", color="blue", label="Antigenically Non-Neutral Mutations")

plot.legend(loc="upper left")

plt.savefig('Antigenic_Fig.png')
