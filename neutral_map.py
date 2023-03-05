import random as rand

AMINO_CHARS = ["A", "C", "G", "U"]

EQUIVALENT_SEQUENCES = [
    ({"GUU", "GUC", "GUA", "GUG"}, "val"),
    ({"GCU", "GCC", "GCA", "GCG"}, "ala"),
    ({"GAU", "GAC"}, "asp"),
    ({"GAA", "GAG"}, "glu"),
    ({"GGU", "GGC", "GGA", "GGG"}, "gly"),
    ({"UUU", "UUC"}, "phe"),
    ({"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"}, "leu"),
    ({"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"}, "ser"),
    ({"UAU", "UAC"}, "tyr"),
    ({"UAA", "UAG", "UGA"}, "stop"),
    ({"UGU", "UGC"}, "cys"),
    ({"UGG"}, "trp"),
    ({"CCU", "CCC", "CCA", "CCG"}, "pro"),
    ({"CAU", "CAC"}, "his"),
    ({"CAA", "CAG"}, "gln"),
    ({"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"}, "arg"),
    ({"AUU", "AUC", "AUA"}, "ile"),
    ({"AUG"}, "met"),
    ({"ACU", "ACC", "ACA", "ACG"}, "thr"),
    ({"AAU", "AAC"}, "asn"),
    ({"AAA", "AAG"}, "lys"),
]


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


def permute_coding_sequence(coding_sequence):
    permutations = []
    for i in range(len(coding_sequence)):
        before = coding_sequence[0:i]
        char = coding_sequence[i]
        after = coding_sequence[i + 1:]
        for new_char in AMINO_CHARS:
            if new_char != char:
                permutations.append(before + new_char + after)
    return permutations


def get_coded_for(sequence):
    for seq in EQUIVALENT_SEQUENCES:
        if sequence in seq[0]:
            return seq[1]


def is_neutral(original_sequence, new_sequence):
    return get_coded_for(original_sequence) == get_coded_for(new_sequence)


def get_all_equivalents(original_sequence, new_sequences):
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
        neutral_seqs = get_all_equivalents(seq, new_seqs)
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
    neutrals = get_all_equivalents("GUU", perms)
    all_new, neutrals = get_new_genomes_and_neutral_genomes("GUUCAAGCA")
    print("krewl")
