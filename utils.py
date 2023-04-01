import pickle as pkl
import math

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


EQUIVALENT_SEQUENCES_ORIG = [
    ({"GUU", "GUC", "GUA", "GUG"}, "V"),
    ({"GCU", "GCC", "GCA", "GCG"}, "A"),
    ({"GAU", "GAC"}, "D"),
    ({"GAA", "GAG"}, "E"),
    ({"GGU", "GGC", "GGA", "GGG"}, "G"),
    ({"UUU", "UUC"}, "F"),
    ({"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"}, "L"),
    ({"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"}, "S"),
    ({"UAU", "UAC"}, "Y"),
    ({"UAA", "UAG", "UGA"}, "*"),
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

EQUIVALENT_SEQUENCES = [(translate_sequence_to_numbers(tup[0]), tup[1]) for tup in EQUIVALENT_SEQUENCES_ORIG]


def translate_numbers_to_string(to_translate):
    ret = ""
    for char in to_translate:
        ret += number_to_char(char)
    return ret


def translate_sequence_to_numbers(to_translate):
    return set([translate_string_to_numbers(trans) for trans in to_translate])


def load_sequence(TRANSLATE_FLAG=True):
    with open("data/spike_protein_genes.txt", "r") as f:
        untranslated = f.read().replace(" ", "").replace("\n", "").strip().upper()
    if TRANSLATE_FLAG:
        return tranlate_proteints_to_base_pairs(untranslated)
    else:
        return untranslated


def tranlate_proteints_to_base_pairs(proteins, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    ret = []
    for protein in proteins:
        found = 0
        for base_pair, prot in translation_list:
            if prot == protein:
                ret.append(list(base_pair)[0])
                found = 1
        if not found:
            raise ValueError(f"PROTEIN {protein} DIDN'T MATCH")
    return "".join(ret)


def translate_base_pairs_to_proteins(base_pairs_to_translate, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    for base_pairs, prot in translation_list:
        if base_pairs_to_translate in base_pairs:
            return prot
    raise ValueError(f"INVALID BASE PAIR SEQUENCE {base_pairs_to_translate} DIDN'T MATCH")

def translate_amino_sequence_to_proteins(amino_sequence_to_translate, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    ret = []
    amino_string = str(amino_sequence_to_translate)
    for i in range(0, len(amino_string), 3):
        local_amino_sequence = amino_string[i:i+3]
        ret.append(translate_base_pairs_to_proteins(local_amino_sequence, translation_list))
    return "".join(ret)


CODED_FOR_DICT = {}


def get_coded_for(sequence, translation_list=None):
    if translation_list is None:
        translation_list = EQUIVALENT_SEQUENCES
    if sequence in CODED_FOR_DICT.keys():
        return CODED_FOR_DICT[sequence]
    for seq in translation_list:
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


NEUTRAL_PROBABILITY_LOW = {
    "V": 1 / 3, "A": 1 / 3, "D": 1 / 9,
    "E": 1 / 9, "G": 1 / 3, "F": 1 / 9,
    "L": 2 / 9, "S": 1 / 9, "Y": 1 / 9,
    "X": 1 / 9, "C": 1 / 9, "W": 0,
    "P": 1 / 3, "H": 2 / 9, "Q": 2 / 9,
    "R": 2 / 9, "I": 2 / 9, "M": 0,
    "T": 1 / 3, "N": 1 / 9, "K": 1 / 9
}

NEUTRAL_PROBABILITY_HIGH = {
    "V": 1 / 3, "A": 1 / 3, "D": 1 / 9,
    "E": 1 / 9, "G": 1 / 3, "F": 1 / 9,
    "L": 4 / 9, "S": 1 / 3, "Y": 1 / 9,
    "X": 2 / 9, "C": 1 / 9, "W": 0,
    "P": 1 / 3, "H": 2 / 9, "Q": 2 / 9,
    "R": 4 / 9, "I": 2 / 9, "M": 0,
    "T": 1 / 3, "N": 1 / 9, "K": 1 / 9
}


def get_neutral_probability(sequence, LOW_FLAG=True):
    probability = 0
    for i in range(0, len(sequence)):
        if LOW_FLAG:
            probability += NEUTRAL_PROBABILITY_LOW[sequence[i]]
        else:
            probability += NEUTRAL_PROBABILITY_HIGH[sequence[i]]
    return probability / len(sequence)


def n_choose_k(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))


def get_approximate_neutral_mutations(sequence, hop):
    low_prob = get_neutral_probability(sequence, True)
    high_prob = get_neutral_probability(sequence, False)
    n = len(sequence) * 3
    k = hop
    total_mutations = n_choose_k(n, k) * 3 ** k
    neutral_mutations_low = total_mutations * low_prob ** k
    neutral_mutations_high = total_mutations * high_prob ** k
    return neutral_mutations_low, neutral_mutations_high, total_mutations


def get_neutral_neighbors(sequence, hop, prob):
    n = len(sequence) * 3
    k = hop

    N = n * 3
    print(N)
    for i in range(1, k + 1):
        N = N - (n_choose_k(n, i) * 3 ** i * prob ** i) + (n_choose_k(n, i + 1) * 3 ** (i + 1) * prob ** i)
        print(N)
    return N


def get_synonymous_nonsynonymous_mutations(sequence, hop):
    low_prob = get_neutral_probability(sequence, True)
    high_prob = get_neutral_probability(sequence, False)

    n = len(sequence) * 3
    k = hop

    N_low = get_neutral_neighbors(sequence, hop, low_prob)
    N_high = get_neutral_neighbors(sequence, hop, high_prob)
    N_synonymous_low = n_choose_k(n, k + 1) * 3 ** (k + 1) * low_prob ** (k + 1)
    N_synonymous_high = n_choose_k(n, k + 1) * 3 ** (k + 1) * high_prob ** (k + 1)
    N_nonsynonymous_low = N_low - N_synonymous_low
    N_nonsynonymous_high = N_high - N_synonymous_high
    return N_synonymous_low, N_nonsynonymous_low, N_synonymous_high, N_nonsynonymous_high


def get_antigenic_synonymous_nonsynonymous_mutations(sequence, hop):
    prob = 0.97409326424

    n = len(sequence) * 3
    k = hop

    N = get_neutral_neighbors(sequence, hop, prob)
    N_synonymous = n_choose_k(n, k + 1) * 3 ** (k + 1) * prob ** (k + 1)
    N_nonsynonymous = N - N_synonymous
    return N_synonymous, N_nonsynonymous

# import pickle as pkl
# with open("results/antigenically_neutral_2_hop_map.pkl", "rb") as f:
#   loaded = pkl.load(f)
#   print(loaded)
