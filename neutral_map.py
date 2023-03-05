import random as rand

AMINO_CHARS = ["A", "C", "G", "U"]

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

def load_sequence():
    with open("data/spike_protein_genes.txt", "r") as f:
        untranslated = f.read().replace(" ","").replace("\n", "").strip().upper()
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

def break_into_coding_seqs(overall_genetic_sequences):
    # looping till length l
    for i in range(0, len(overall_genetic_sequences), 3):
        yield overall_genetic_sequences[i:i + 3]


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
    for i in range(0, len(genetic_sequence), 3):
        before = genetic_sequence[0:i]
        seq = genetic_sequence[i: i + 3]
        after = genetic_sequence[i + 3:]
        new_seqs = permute_coding_sequence(seq)
        neutral_seqs = extract_all_equivalents(seq, new_seqs)
        this_round = [before + new_seq + after for new_seq in new_seqs]
        this_round_neutral = [before + neutral_seq + after for neutral_seq in neutral_seqs]
        new_sequences.extend(this_round)
        neutral_sequences.extend(this_round_neutral)
    return new_sequences, neutral_sequences


TEST_SEQUENCE = "GUUCAAGCA"

if __name__ == "__main__":
    print(generate_genome(9))
    check_equivalents()
    perms = permute_coding_sequence("GUU")
    equivalents = extract_all_equivalents("GUU", perms)
    all_new, neutrals = get_new_genomes_and_neutral_genomes("GUUCAAGCA")
    l = load_sequence()
    all_new_loaded, neutrals_loaded = get_new_genomes_and_neutral_genomes(l)
    results = [set(neutrals_loaded)]
    for i in range(1,5):
        loop_neuts = set()
        j = 0
        for sequence in results[i - 1]:
            all_new_loaded_loop, neutrals_loaded_loop = get_new_genomes_and_neutral_genomes(sequence)
            new = set(neutrals_loaded_loop) - results[i - 1]
            loop_neuts.update(neutrals_loaded_loop)
            j += 1
            if j % 100 == 0:
                print(j)
        results.append(loop_neuts)
        print("krewl")
