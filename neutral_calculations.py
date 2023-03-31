import matplotlib.pyplot as plt
from matplotlib import style
from utils import *

sequence = load_sequence(False)
synonymous_low = []
synonymous_high = []
nonsynonymous_low = []
nonsynonymous_high = []
print("Probability of a single mutation away from original genome resulting in a neutral mutation: \n"+str(round(get_neutral_probability(sequence) * 100, 3))+"% <= probability <= "+str(round(get_neutral_probability(sequence, False) * 100, 3))+"%\n")
xs = list(range(1, 11))
for hop in xs:
    syn_low, nonsyn_low, syn_high, nonsyn_high = get_synonymous_nonsynonymous_mutations(sequence, hop)
    synonymous_low.append(syn_low)
    synonymous_high.append(syn_high)
    nonsynonymous_low.append(nonsyn_low)
    nonsynonymous_high.append(nonsyn_high)

style.use('fivethirtyeight')

plt.rc('axes', titlesize=14)
plt.rc('axes', labelsize=12)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)

fig = plt.figure(figsize=(8, 6))
plot = fig.add_subplot(1, 1, 1)

plot.set_title("Neutral Mutations vs Non-Neutral Mutations in Neighborhood at n-hop")

plot.set_xlabel("n (steps)")
plot.set_xticks(range(0, int(max(xs))+1, 1))

plot.set_ylabel("unique mutations in neighborhood")
plot.set_yscale("log")

plt.tight_layout()

syn_err = [(y-x)/2 for x,y in zip(synonymous_low, synonymous_high)]
nonsyn_err = [(y-x)/2 for x,y in zip(nonsynonymous_low, nonsynonymous_high)]
average_syn = [(x+y)/2 for x,y in zip(synonymous_low, synonymous_high)]
average_nonsyn = [(x+y)/2 for x,y in zip(nonsynonymous_low, nonsynonymous_high)]

plot.fill_between(xs, synonymous_low, synonymous_high, alpha=0.2, color="red")
#plot.errorbar(xs, average_neutral, yerr=err, marker="o", capsize=4, linestyle="dashed", linewidth=1, color="red", label="Neutral Mutations")
plot.plot(xs, average_syn, linewidth=1, marker="o", color="red", label="Neutral Mutations")

plot.fill_between(xs, nonsynonymous_low, nonsynonymous_high, alpha=0.2, color="blue")
#plot.errorbar(xs, average_neutral, yerr=err, marker="o", capsize=4, linestyle="dashed", linewidth=1, color="red", label="Neutral Mutations")
plot.plot(xs, average_nonsyn, linewidth=1, marker="o", color="blue", label="Non-Neutral Mutations")

plot.legend(loc="upper left")

plt.savefig('Neutral_Fig.png')
