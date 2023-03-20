import matplotlib.pyplot as plt
from matplotlib import style
from utils import *

sequence = load_sequence(False)
total = []
low_prob = []
high_prob = []
print("Probability of a single mutation away from original genome resulting in a neutral mutation: \n"+str(round(get_neutral_probability(sequence) * 100, 3))+"% <= probability <= "+str(round(get_neutral_probability(sequence, False) * 100, 3))+"%\n")
for hop in range(0, 6):
    low, high, tot = get_approximate_neutral_mutations(sequence, hop)
    total.append(tot)
    low_prob.append(low)
    high_prob.append(high)

xs = [0, 1, 2, 3, 4, 5]

style.use('fivethirtyeight')

plt.rc('axes', titlesize=12)
plt.rc('axes', labelsize=10)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)

fig = plt.figure(figsize=(8, 6))
plot = fig.add_subplot(1, 1, 1)

plot.set_title("Computed Neutral Mutations vs Total Mutations of Spike Protein")

plot.set_xlabel("n (steps)")
plot.set_xticks(range(0, int(max(xs))+1, 1))

plot.set_ylabel("unique mutations")
plot.set_yscale("log")

plt.tight_layout()

plot.plot(xs, total, marker="o", linestyle="dashed", linewidth=1, color="blue", label="Total Mutations")

#plot.fill_between(xs, low_prob, high_prob, alpha=0.5, color="red")
#plot.plot(xs, [(x+y)/2 for x,y in zip(low_prob, high_prob)], marker="o", linestyle="dashed", linewidth=1, color="red", label="Neutral Mutations")

plot.fill_between(xs, low_prob, high_prob, alpha=0.2, color="red")
plot.errorbar(xs, [(x+y)/2 for x,y in zip(low_prob, high_prob)], yerr=[(y-x)/2 for x,y in zip(low_prob, high_prob)], marker="o", capsize=4, linestyle="dashed", linewidth=1, color="red", label="Neutral Mutations")

plot.legend(loc="upper left")

plt.savefig('Neutral_Fig.png')
