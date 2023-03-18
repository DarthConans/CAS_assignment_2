import multiprocessing
import random as rand
import pickle as pkl
import os
from multiprocessing import Pool
from utils import *
from bindingcalculator import *




def neutral_genomes_and_sites(genetic_sequence):
    neutral_sequences = []
    translated_sequence = str(genetic_sequence)
    for i in range(0, len(translated_sequence), 3):
        before = translated_sequence[0:i]
        seq = translated_sequence[i: i + 3]
        after = translated_sequence[i + 3:]
        new_seqs = permute_coding_sequence(seq)
        neutral_seqs = extract_all_equivalents(seq, new_seqs)
        this_round_neutral = [(int(before + neutral_seq + after), i // 3) for neutral_seq in neutral_seqs]
        neutral_sequences.extend(this_round_neutral)
    return neutral_sequences

ESCAPE_THRESHOLD = .1

def antigenic_neutral_function(antigen_frame):
    original = antigen_frame["original_escape"].sum()
    retained = antigen_frame["retained_escape"].sum()
    difference = abs(original - retained)
    if retained > original or difference < ESCAPE_THRESHOLD:
        return True
    return False

binding = BindingCalculator()



def get_antigenically_neutral(neutral_tuples, offset = 331):
    ret = []
    sites_have = set([neutral_tuple[1] + offset for neutral_tuple in neutral_tuples])
    sites_wanted = binding.sites
    count_hits = len(sites_wanted.intersection(sites_have))
    expected = len(neutral_tuples)
    done = 0
    for neutral_tuple in neutral_tuples:
        if len(neutral_tuple) == 2:
            sites = {neutral_tuple[1] + offset}
        else:
            sites = {neutral_tuple[1] + offset, neutral_tuple[2] + offset}
        if len(sites - binding.sites) == 0:
            antigentic_frame = binding.escape_per_site(sites)
            if antigenic_neutral_function(antigentic_frame):
                ret.append(neutral_tuple)
        done += 1
        if done % 1000 == 0:
            print(f"FINISHED {done} OUT OF {expected}")
    return [rets[0] for rets in ret]



if __name__ == "__main__":
    # all_new, neutrals = get_new_genomes_and_neutral_genomes(translate_string_to_numbers("GUUCAAGCA"))
    l = int(load_sequence())
    assert len(str(l)) == 918
    neutrals_loaded = neutral_genomes_and_sites(l)
    neutrals = set(neutrals_loaded)
    antigenically_neutral_1_hop_map = get_antigenically_neutral(neutrals)
    with open("results/antigenically_neutral_1_hop_map.pkl", "wb") as f:
        pkl.dump(antigenically_neutral_1_hop_map, f)
    neutral_2_hop = set()
    for l in neutrals:
        to_neutrify = l[0]
        r = neutral_genomes_and_sites(to_neutrify)
        r = [(a[0], a[1], l[1]) for a in r]
        neutral_2_hop.update(r)

    antigenically_neutral_2_hop_map = get_antigenically_neutral(neutral_2_hop)


    with open("results/antigenically_neutral_2_hop_map.pkl", "wb") as f:
        pkl.dump(antigenically_neutral_2_hop_map, f)
    # with open("results/loop_neuts_hop_3.pkl", "rb") as f:

    #    neutrals_loaded = pkl.load(f)
    # neutrals_loaded = list(break_into_chunks(list(neutrals_loaded), 100000))
    # i = 0
    # total = len(neutrals_loaded)
    # for neutral_chunk in neutrals_loaded:
    #    with open(f"broken_up/hop_3_{i}.pkl", "wb") as f:
    #        pkl.dump(neutral_chunk, f)
    #    i += 1
    #    print(f"I'VE DUMPED {i} OF {total}")
    # results = [set(neutrals_loaded)]
    # with open(f"results/loop_neuts_hop_1.pkl", "wb") as f:
    #    pkl.dump(set(neutrals_loaded), f)
    # for i in range(1,5):
    # loop_neuts = set()
    files = os.listdir("broken_up")
    sorted(files)
    print("krewl")
