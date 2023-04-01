import multiprocessing
import random as rand
import pickle as pkl
import os
from multiprocessing import Pool
from utils import *







if __name__ == "__main__":
    #generate_antigenic_from_sites()
    l = int(load_sequence())
    assert len(str(l)) == 579
    neutrals_loaded = neutral_genomes_and_sites(l)
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
        print(f"I'VE DONE {done} OUT OF {to_do}")


    antigenically_neutral_2_hop_map = get_antigenically_neutral(neutral_2_hop_all)


    with open("results/antigenically_neutral_2_hop_all_map.pkl", "wb") as f:
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
