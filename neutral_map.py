import multiprocessing
import random as rand
import pickle as pkl
import os
from utils import *
from multiprocessing import Pool


TEST_SEQUENCE = "GUUCAAGCA"


def process_broken_up_file(file_name):
    with open(f"broken_up/{file_name}", "rb") as f:
        sequences = pkl.load(f)
    file_neuts = set()
    j = 0
    total = len(sequences)
    for sequence in sequences:
        all_new_loaded_loop, neutrals_loaded_loop = get_new_genomes_and_neutral_genomes(sequence)

        file_neuts.update(neutrals_loaded_loop)
        j += 1
        if j % 100 == 0:
            print(f"{j} OUT OF {total} DONE ON FILE {file_name}. I HAVE {len(file_neuts)} RESULTS")
    with open(f"results/hop4/{file_name}", "wb") as f:
        pkl.dump(file_neuts, f)


def serialize_results(neutral_mutations, hop):
    with open(f"results/loop_neuts_hop_{hop}.pkl", "wb") as f:
        pkl.dump(neutral_mutations, f)


if __name__ == "__main__":
    # all_new, neutrals = get_new_genomes_and_neutral_genomes(translate_string_to_numbers("GUUCAAGCA"))
    l = int(load_sequence())
    assert len(str(l)) == 579
    _, neutrals_loaded = get_new_genomes_and_neutral_genomes(l)
    neutrals = set(neutrals_loaded)
    serialize_results(neutrals_loaded, 1)
    last_loop_neuts = neutrals
    overall_neuts = neutrals
    for i in range(2, 4):
        this_loop_neuts = set()
        done = 0
        to_do = len(last_loop_neuts)
        for l in last_loop_neuts:
            result = get_new_genomes_and_neutral_genomes(l)
            this_loop_neuts.update(result[1])
            done += 1
            if done % 100 == 0:
                print(f"I'VE DONE {done} OF {to_do} ON HOP {i}")
        this_loop_neuts = this_loop_neuts - overall_neuts
        overall_neuts.update(this_loop_neuts)
        serialize_results(this_loop_neuts, i)
        last_loop_neuts = this_loop_neuts.copy()
        print("krewl")
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
