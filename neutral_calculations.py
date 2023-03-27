import matplotlib.pyplot as plt
from matplotlib import style
from utils import *

sequence = load_sequence(False)
total = []
low_prob = []
high_prob = []
print("Probability of a single mutation away from original genome resulting in a neutral mutation: \n"+str(round(get_neutral_probability(sequence) * 100, 3))+"% <= probability <= "+str(round(get_neutral_probability(sequence, False) * 100, 3))+"%\n")
xs = list(range(1, 5))
for hop in xs:
    low, high, tot = get_approximate_neutral_mutations(sequence, hop)
    total.append(tot)
    low_prob.append(low)
    high_prob.append(high)

style.use('fivethirtyeight')

plt.rc('axes', titlesize=14)
plt.rc('axes', labelsize=12)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)

fig = plt.figure(figsize=(8, 6))
plot = fig.add_subplot(1, 1, 1)

plot.set_title("Computed Neutral Mutations vs Total Mutations of Spike Protein")

plot.set_xlabel("n (steps)")
plot.set_xticks(range(0, int(max(xs))+1, 1))

plot.set_ylabel("unique mutations")
plot.set_yscale("log")

plt.tight_layout()

plot.plot(xs, total, marker="o", linestyle="dashed", linewidth=1, color="purple", label="Total Mutations")
#plot.plot(xs, total, linewidth=1, color="purple", label="Total Mutations")

#plot.fill_between(xs, low_prob, high_prob, alpha=0.5, color="red")
#plot.plot(xs, [(x+y)/2 for x,y in zip(low_prob, high_prob)], marker="o", linestyle="dashed", linewidth=1, color="red", label="Neutral Mutations")

err = [(y-x)/2 for x,y in zip(low_prob, high_prob)]
average_neutral = [(x+y)/2 for x,y in zip(low_prob, high_prob)]
average_non_neutral = [x-y for x,y in zip(total, average_neutral)]

plot.fill_between(xs, low_prob, high_prob, alpha=0.2, color="red")
#plot.errorbar(xs, average_neutral, yerr=err, marker="o", capsize=4, linestyle="dashed", linewidth=1, color="red", label="Neutral Mutations")
plot.plot(xs, average_neutral, linewidth=1, marker="o", color="red", label="Neutral Mutations")

plot.fill_between(xs, [x-y for x,y in zip(average_non_neutral, err)], [x+y for x,y in zip(average_non_neutral, err)], alpha=0.2, color="blue")
#plot.errorbar(xs, average_non_neutral, yerr=err, marker="o", capsize=4, linestyle="dashed", linewidth=1, color="blue", label="Non-Neutral Mutations")
plot.plot(xs, average_non_neutral, linewidth=1, marker="o", color="blue", label="Non-Neutral Mutations")

plot.legend(loc="upper left")

plt.savefig('Neutral_Fig.png')
