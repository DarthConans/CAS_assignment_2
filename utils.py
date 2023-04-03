import pickle as pkl
import math

from bindingcalculator import *

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
    for char in str(to_translate):
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

def all_genomes_and_sites(genetic_sequence):
    sequences = set()
    translated_sequence = str(genetic_sequence)
    for i in range(0, len(translated_sequence), 3):
        before = translated_sequence[0:i]
        seq = translated_sequence[i: i + 3]
        after = translated_sequence[i + 3:]
        new_seqs = permute_coding_sequence(seq)
        this_round_new = [(int(before + seq + after), i // 3) for seq in new_seqs]
        sequences.update(this_round_new)
    return sequences

ESCAPE_THRESHOLD = 3

def antigenic_neutral_function(antigen_frame):
    original = antigen_frame["original_escape"].sum()
    retained = antigen_frame["retained_escape"].sum()
    difference = abs(original - retained)
    if retained > original or difference < ESCAPE_THRESHOLD:
        return True
    return False

def get_antigen_value(antigen_frame):
    original = antigen_frame["original_escape"].sum()
    retained = antigen_frame["retained_escape"].sum()
    difference = abs(original - retained)
    return original, retained, difference


binding = BindingCalculator()

non_antigenic = set()

anitgen_memo = {

}
def get_antigenically_neutral(neutral_tuples, offset=330):
    ret = []
    sites_have = set([neutral_tuple[1] + offset for neutral_tuple in neutral_tuples])
    sites_wanted = binding.sites
    count_hits = len(sites_wanted.intersection(sites_have))
    expected = len(neutral_tuples)
    done = 0
    for neutral_tuple in neutral_tuples:
        if len(neutral_tuple) == 2:
            sites = (neutral_tuple[1] + offset,)
        else:
            sites = tuple(sorted((neutral_tuple[1] + offset, neutral_tuple[2] + offset)))
        if sites in anitgen_memo.keys():
            truth = anitgen_memo[sites]
        else:

            truth = False
            site_set = set()
            for site in sites:
                site_set.add(site)
            if len(site_set - binding.sites) == 0:
                antigentic_frame = binding.escape_per_site(sites)
                truth = antigenic_neutral_function(antigentic_frame)
            anitgen_memo[sites] = truth
        if truth:
            ret.append(neutral_tuple)
        done += 1
        if done % 1000 == 0:
            print(f"FINISHED {done} OUT OF {expected}")
    return [(rets[0], rets[1]) for rets in ret]

def get_antigenically_non_neutral(neutral_tuples, offset=330):
    ret = []
    sites_have = set([neutral_tuple[1] + offset for neutral_tuple in neutral_tuples])
    sites_wanted = binding.sites
    count_hits = len(sites_wanted.intersection(sites_have))
    expected = len(neutral_tuples)
    done = 0
    antigenically_active = {345, 443, 483, 485, 489}
    for neutral_tuple in neutral_tuples:
        if len(neutral_tuple) == 2:
            sites = {neutral_tuple[1] + offset}
        else:
            sites = {neutral_tuple[1] + offset, neutral_tuple[2] + offset}
        t = any([site for site in sites if site in antigenically_active])
        if len(sites - binding.sites) == 0 and any([site for site in sites if site in antigenically_active]):

            antigentic_frame = binding.escape_per_site(sites)
            if not antigenic_neutral_function(antigentic_frame):
                values = get_antigen_value(antigentic_frame)
                ret.append((neutral_tuple, (values[0], values[1], values[2])))
        done += 1
        if done % 1000 == 0:
            print(f"FINISHED {str(done)} OUT OF {str(expected)}")
            print(f"RET IS CURRENTLY {len(ret)} LONG")
    return ret


def generate_antigenic_from_sites():
    one_hop = set()
    for site in binding.sites:
        frame = binding.escape_per_site([site])
        if antigenic_neutral_function(frame):
            one_hop.add(site)
    with open("results/only_antigenically_neutral_1_hop_map.pkl", "wb") as f:
        pkl.dump(one_hop, f)
    two_hops = set()
    to_do = len(one_hop)
    done = 0
    for first_site in one_hop:
        for second_site in binding.sites:
            if first_site != second_site:
                frame = binding.escape_per_site([first_site, second_site])
                if antigenic_neutral_function(frame):
                    two_hops.add((first_site, second_site))
        done += 1
        print(f"I'VE DONE {done} OUT OF {to_do}")
    with open("results/only_antigenically_neutral_2_hop_map.pkl", "wb") as f:
        pkl.dump(one_hop, f)
    return one_hop, two_hops


def probability_epistatic_compensatory_mutation(sequence, hop):
    p2 = 0.00001
    p_low = get_neutral_probability(sequence, True)
    p_high = get_neutral_probability(sequence, False)

    n = len(sequence) * 3
    k = hop

    # E_0
    E_0_low = n_choose_k(n, 1) * 3 * p2 * (1-p_low)
    E_0_high = n_choose_k(n, 1) * 3 * p2 * (1-p_high)
    E_low = 0
    E_high = 0
    for i in range(1, k+1):
        #E_n
        E_low = E_low + (n_choose_k(n, i+1) * (3**(i+1)) * (p_low**i) * (1-p_low) * ((p2*((p2**i)-((1-p2)**i)))/(2*p2 - 1)))
        E_high = E_high + (n_choose_k(n, i+1) * (3**(i+1)) * (p_high**i) * (1-p_high) * ((p2*((p2**i)-((1-p2)**i)))/(2*p2 - 1)))
    E_low = E_0_low + E_low
    E_high = E_0_high + E_high

    S_n_low, S_n_bar_low, S_n_high, S_n_bar_high = get_synonymous_nonsynonymous_mutations(sequence, hop)

    print('En low', E_low, S_n_low, S_n_bar_low)
    print('En high', E_high, S_n_high, S_n_bar_high)

    Pe_low = E_low / S_n_bar_low
    Pe_high = E_high / S_n_bar_high
    return Pe_low, Pe_high


# import pickle as pkl
# with open("results/antigenically_neutral_2_hop_map.pkl", "rb") as f:
#   loaded = pkl.load(f)
#   print(loaded)
