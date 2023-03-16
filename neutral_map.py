import multiprocessing
import random as rand
import pickle as pkl
import os
from multiprocessing import Pool

AMINO_CHARS = ["1", "2", "3", "4"]

translate_dict = {
    "G": "1",
    "U": "2",
    "C": "3",
    "A": "4",
}


def char_to_number(char):
    return translate_dict[char]


translate_back_dict = {
    "1": "G",
    "2": "U",
    "3": "C",
    "4": "A",
}


def number_to_char(char):
    return translate_back_dict[char]


def translate_string_to_numbers(to_translate):
    ret = ""
    for char in to_translate:
        ret += char_to_number(char)
    return ret


def translate_sequence_to_numbers(to_translate):
    return set([translate_string_to_numbers(trans) for trans in to_translate])


EQUIVALENT_SEQUENCES = [
    ({"GUU", "GUC", "GUA", "GUG"}, "V"),
    ({"GCU", "GCC", "GCA", "GCG"}, "A"),
    ({"GAU", "GAC"}, "D"),
    ({"GAA", "GAG"}, "E"),
    ({"GGU", "GGC", "GGA", "GGG"}, "G"),
    ({"UUU", "UUC"}, "F"),
    ({"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"}, "L"),
    ({"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"}, "S"),
    ({"UAU", "UAC"}, "Y"),
    ({"UAA", "UAG", "UGA"}, "X"),
    ({"UGU", "UGC"}, "C"),
    ({"UGG"}, "W"),
    ({"CCU", "CCC", "CCA", "CCG"}, "P"),
    ({"CAU", "CAC"}, "H"),
    ({"CAA", "CAG"}, "Q"),
    ({"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"}, "R"),
    ({"AUU", "AUC", "AUA"}, "I"),
    ({"AUG"}, "M"),
    ({"ACU", "ACC", "ACA", "ACG"}, "T"),
    ({"AAU", "AAC"}, "N"),
    ({"AAA", "AAG"}, "K"),
]

EQUIVALENT_SEQUENCES = [(translate_sequence_to_numbers(tup[0]), tup[1]) for tup in EQUIVALENT_SEQUENCES]


def translate_numbers_to_string(to_translate):
    ret = ""
    for char in to_translate:
        ret += number_to_char(char)
    return ret


def load_sequence():
    with open("data/spike_protein_genes.txt", "r") as f:
        untranslated = f.read().replace(" ", "").replace("\n", "").strip().upper()
    return tranlate_proteints_to_base_pairs(untranslated)


def tranlate_proteints_to_base_pairs(proteins):
    ret = []
    for protein in proteins:
        found = 0
        for base_pair, prot in EQUIVALENT_SEQUENCES:
            if prot == protein:
                ret.append(list(base_pair)[0])
                found = 1
        if not found:
            raise ValueError(f"PROTEIN {protein} DIDN'T MATCH")
    return "".join(ret)


def break_into_chunks(to_chunkify, count):
    # looping till length l
    for i in range(0, len(to_chunkify), count):
        yield to_chunkify[i:i + count]


def generate_genome(length=6):
    if length < 0 or length % 3 != 0:
        raise ValueError(f"{length} WON'T WORK. IT NEEDS TO BE GREATER THAN 0 AND DIVISIBLE BY 3.")
    acids = rand.choices(AMINO_CHARS, k=length)
    return "".join(acids)


def check_equivalents():
    coded_for = [s[1] for s in EQUIVALENT_SEQUENCES]
    assert len(coded_for) == len(set(coded_for))
    amino_combos = []
    for amino_combo in EQUIVALENT_SEQUENCES:
        amino_combos.extend(amino_combo[0])
    assert len(amino_combos) == len(set(amino_combos)) == (4 * 4 * 4)


PERMUTATION_DICT = {}


def permute_coding_sequence(coding_sequence):
    if coding_sequence in PERMUTATION_DICT.keys():
        return PERMUTATION_DICT[coding_sequence]
    permutations = []
    for i in range(len(coding_sequence)):
        before = coding_sequence[0:i]
        char = coding_sequence[i]
        after = coding_sequence[i + 1:]
        for new_char in AMINO_CHARS:
            if new_char != char:
                permutations.append(before + new_char + after)
    PERMUTATION_DICT[coding_sequence] = permutations
    return permutations


CODED_FOR_DICT = {}


def get_coded_for(sequence):
    if sequence in CODED_FOR_DICT.keys():
        return CODED_FOR_DICT[sequence]
    for seq in EQUIVALENT_SEQUENCES:
        if sequence in seq[0]:
            CODED_FOR_DICT[sequence] = seq[1]
            return seq[1]


def is_neutral(original_sequence, new_sequence):
    return get_coded_for(original_sequence) == get_coded_for(new_sequence)


def extract_all_equivalents(original_sequence, new_sequences):
    original_coded_for = get_coded_for(original_sequence)
    return [seq for seq in new_sequences if get_coded_for(seq) == original_coded_for]


def get_new_genomes_and_neutral_genomes(genetic_sequence):
    new_sequences = []
    neutral_sequences = []
    translated_sequence = str(genetic_sequence)
    for i in range(0, len(translated_sequence), 3):
        before = translated_sequence[0:i]
        seq = translated_sequence[i: i + 3]
        after = translated_sequence[i + 3:]
        new_seqs = permute_coding_sequence(seq)
        neutral_seqs = extract_all_equivalents(seq, new_seqs)
        this_round = [int(before + new_seq + after) for new_seq in new_seqs]
        this_round_neutral = [int(before + neutral_seq + after) for neutral_seq in neutral_seqs]
        new_sequences.extend(this_round)
        neutral_sequences.extend(this_round_neutral)
    return new_sequences, neutral_sequences


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
    assert len(str(l)) == 918
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
    #with open("results/loop_neuts_hop_3.pkl", "rb") as f:

    #    neutrals_loaded = pkl.load(f)
    #neutrals_loaded = list(break_into_chunks(list(neutrals_loaded), 100000))
    #i = 0
    #total = len(neutrals_loaded)
    #for neutral_chunk in neutrals_loaded:
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
