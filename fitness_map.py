import os.path

import numpy as np
import pandas as pd
import pickle as pkl
from utils import *


def get_mutations_and_indexes(original_genome, new_genome):
    ret = frozenset([(i, original_genome[i], new_genome[i]) for i in range(len(original_genome)) if original_genome[i] != new_genome[i]])
    return ret

FITNESS_FRAME = pd.read_csv("data/aamut_fitness_all.csv")[["clade_founder_aa", "mutant_aa", "aa_site", "delta_fitness"]]

def get_fitness_change(mutation_data, offset=330):
    total_change = 0
    for mutation_datum in mutation_data:
        site = mutation_datum[0] + offset
        original_amino_acid = mutation_datum[1]
        new_amino_acid = mutation_datum[2]
        fitness_row = FITNESS_FRAME[(FITNESS_FRAME["aa_site"] == site) & (FITNESS_FRAME["clade_founder_aa"] == original_amino_acid) & (FITNESS_FRAME["mutant_aa"] == new_amino_acid)]
        if not fitness_row.empty:
            total_change += fitness_row["delta_fitness"].iloc[0]
        else:
            return np.NAN
    return total_change

def to_csv(path, data, offset=330):
    header = "site,old amino,new amino,fitness change\n"
    row_data = []
    for datum in data:
        row_vals = []
        for frozen in datum[0]:
            site = frozen[0] + offset
            old_amino = frozen[1]
            new_amino = frozen[2]
            row_vals.extend([str(site), old_amino, new_amino])
        while len(row_vals) < 3:
            row_vals.append("None")
        fitness_change = str(datum[1])
        row_vals.append(fitness_change)
        row = ",".join(row_vals)
        row_data.append(row)
    with open(path, "w") as f:
        f.write(header)
        for row in row_data:
            f.write(row + "\n")


if __name__ == "__main__":
    original = load_sequence(False)
    l = int(load_sequence())

    one_hop_all, one_hop_neut = get_new_genomes_and_neutral_genomes(l)
    one_hop_all = set(one_hop_all)
    one_hop_neut = set(one_hop_neut)
    one_hop_all_translated = set([translate_amino_sequence_to_proteins(x) for x in one_hop_all])
    two_hop_all = set()
    for neut in one_hop_neut:
        all_mutations, _ = get_new_genomes_and_neutral_genomes(neut)
        two_hop_all.update(all_mutations)
    if os.path.exists("results/two_hop_translated.pkl"):
        with open("results/two_hop_translated.pkl", "rb") as f:
            two_hop_all_translated = pkl.load(f)
    else:
        two_hop_all_translated = set()
        done = 0
        to_do = len(two_hop_all)
        for sequence in two_hop_all:
            translated = translate_amino_sequence_to_proteins(sequence)
            two_hop_all_translated.add(translated)
            done += 1
            if done % 1000 == 0:
                print(f"I'VE TRANSLATED {done} OUT OF {to_do}")
                print(f"THERE ARE {len(two_hop_all_translated)} IN THE SET")
        with open("results/two_hop_translated.pkl", "wb") as f:
            pkl.dump(two_hop_all_translated, f)
    one_and_two_all_translated = one_hop_all_translated.union(two_hop_all_translated)
    mutations_and_indexes = set([get_mutations_and_indexes(original, x) for x in two_hop_all_translated])
    mutations_and_indexes = set([x for x in mutations_and_indexes if x])
    mutations_indexes_and_fitness_changes = [(x, get_fitness_change(x)) for x in mutations_and_indexes]
    mutations_indexes_and_fitness_changes = [x for x in mutations_indexes_and_fitness_changes if x[1] == x[1]]
    mutations_indexes_and_fitness_changes = sorted(mutations_indexes_and_fitness_changes, key=lambda x: x[1], reverse=True)
    to_csv("results/fitness_change_1_and_2_hop.csv", mutations_indexes_and_fitness_changes)
    print("krewl")