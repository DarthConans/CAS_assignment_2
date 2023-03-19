from utils import *

sequence = load_sequence(False)
print("Probability of a single mutation away from original genome resulting in a neutral mutation: \n"+str(round(get_neutral_probability(sequence) * 100, 3))+"%\n")
for hop in range(0, 5):
    get_approximate_neutral_mutations(sequence, hop)