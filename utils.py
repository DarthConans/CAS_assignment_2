import pickle as pkl

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


def translate_sequence_to_numbers(to_translate):
    return set([translate_string_to_numbers(trans) for trans in to_translate])


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


CODED_FOR_DICT = {}


def get_coded_for(sequence):
    if sequence in CODED_FOR_DICT.keys():
        return CODED_FOR_DICT[sequence]
    for seq in EQUIVALENT_SEQUENCES:
        if sequence in seq[0]:
            CODED_FOR_DICT[sequence] = seq[1]
            return seq[1]


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


def extract_all_equivalents(original_sequence, new_sequences):
    original_coded_for = get_coded_for(original_sequence)
    return [seq for seq in new_sequences if get_coded_for(seq) == original_coded_for]


def is_neutral(original_sequence, new_sequence):
    return get_coded_for(original_sequence) == get_coded_for(new_sequence)


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


def serialize_results(neutral_mutations, hop):
    with open(f"results/loop_neuts_hop_{hop}.pkl", "wb") as f:
        pkl.dump(neutral_mutations, f)
