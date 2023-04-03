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
    for mutation_datum in mutation_data[0]:
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
    header = "site,old amino,new amino,mutation 1 site,mutation 2 site,old antigen escape value,new antigen escape value,antigen difference,fitness change,old sequence,new sequence\n"
    row_data = []
    for datum in data:
        row_vals = []
        for frozen in datum[0][0]:
            site = frozen[0] + offset
            old_amino = frozen[1]
            new_amino = frozen[2]
            row_vals.extend([str(site), old_amino, new_amino])
        while len(row_vals) < 3:
            row_vals.append("None")
        mutation_1_site = str(datum[0][1] + offset)
        mutation_2_site = str(datum[0][2] + offset)
        old_escape_val = str(datum[0][3])
        new_escape_val = str(datum[0][4])
        difference = str(datum[0][5])
        new_sequence = str(datum[0][6])
        old_sequence = str(datum[0][7])
        fitness_change = str(datum[1])
        row_vals.extend([mutation_1_site, mutation_2_site, old_escape_val, new_escape_val, difference, fitness_change, old_sequence, new_sequence])
        row = ",".join(row_vals)
        row_data.append(row)
    with open(path, "w") as f:
        f.write(header)
        for row in row_data:
            f.write(row + "\n")


if __name__ == "__main__":
    original = load_sequence(TRANSLATE_FLAG=False)
    initial = int(load_sequence())
    assert len(str(initial)) == 579
    neutrals_loaded = neutral_genomes_and_sites(initial)
    neutrals = set(neutrals_loaded)
    neutral_2_hop_all = set()
    done = 0
    to_do = len(neutrals)
    neutrals_seq_only = set([x[0] for x in neutrals])
    for l in neutrals:
        to_neutrify = l[0]
        r = all_genomes_and_sites(to_neutrify)
        this_round_unique = set([x[0] for x in r])
        r = [(a[0], a[1], l[1]) for a in r if a[0] not in neutrals_seq_only]
        neutral_2_hop_all.update(r)
        done += 1
    antigenically_non_neutral_2_hop_map = get_antigenically_non_neutral(list(neutral_2_hop_all))
    two_hop_all_translated = set()
    for entry in antigenically_non_neutral_2_hop_map:
        two_hop_all_translated.add((translate_codon_sequence_to_aas(entry[0][0]), entry[0][1], entry[0][2], entry[1][0], entry[1][1], entry[1][2], entry[0][0], initial))
    mutations_and_indexes = set([(get_mutations_and_indexes(original, x[0]), x[1], x[2], x[3], x[4], x[5], translate_numbers_to_string(x[6]), translate_numbers_to_string(x[7])) for x in two_hop_all_translated])
    mutations_and_indexes = set([x for x in mutations_and_indexes if x[0]])
    mutations_indexes_and_fitness_changes = [(x, get_fitness_change(x)) for x in mutations_and_indexes]
    mutations_indexes_and_fitness_changes = [x for x in mutations_indexes_and_fitness_changes if x[1] == x[1]]
    mutations_indexes_and_fitness_changes = sorted(mutations_indexes_and_fitness_changes, key=lambda x: x[1], reverse=True)
    to_csv("results/fitness_change_of_antigenically_non_neutral_2_hop.csv", mutations_indexes_and_fitness_changes)
    print("krewl")