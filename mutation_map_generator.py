import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as cls
import matplotlib.patches as patches
from utils import *

# create a graph where each node represents
# a unique genetic sequence, and edges
# represent a single mutation difference
def create_graph(sequences, fitness_change, antigen_difference, original_length):
    # create graph
    G = nx.Graph()

    # calculate range
    f_min = min(fitness_change[:original_length])
    f_max = max(fitness_change[:original_length])
    f_range = f_max - f_min
    a_min = min(antigen_difference[:original_length])
    a_max = max(antigen_difference[:original_length])
    a_range = a_max - a_min

    # add nodes
    color_map_fitness = []
    color_map_antigen = []
    for i in range(len(sequences)):
        saturation_f = ((fitness_change[i] - f_min) / f_range) * (0.8) + 0.2
        saturation_a = ((antigen_difference[i] - a_min) / a_range) * (0.8) + 0.2
        value = 0.8
        if antigen_difference[i]==0:
            value = 0
        color_map_fitness.append(cls.to_hex(cls.hsv_to_rgb([0.37, saturation_f, value])))
        color_map_antigen.append(cls.to_hex(cls.hsv_to_rgb([0.0, saturation_a, value])))

        G.add_node(i)

    # add edges
    for i in range(len(sequences)):
        for j in range(i, len(sequences)):
            if i!=j and one_mutation_away(sequences[i], sequences[j]):
                G.add_edge(i, j)
    
    fig = plt.figure(dpi=240)
    fig.add_subplot(1, 2, 1)
    nx.draw(G, pos=nx.spring_layout(G, seed=31), node_size=10, node_color=color_map_fitness, edge_color="gray")
    plt.title("Viral Fitness")
    fig.add_subplot(1, 2, 2)
    nx.draw(G, pos=nx.spring_layout(G, seed=31), node_size=10, node_color=color_map_antigen, edge_color="gray")
    plt.title("Antibody Escape")

    synonymous = patches.Patch(color="black", label="Synonymous Mutations")
    viral_fitness = patches.Patch(color="lightGreen", label="Degree of Viral Fitness")
    antigenic_fitness = patches.Patch(color="red", label="Degree of Antibody Escape")

    plt.suptitle("Connected Network to Selected Nodes")
    plt.rc('legend', fontsize=10)
    plt.legend(handles=[synonymous, viral_fitness, antigenic_fitness], loc="lower left", bbox_to_anchor=(-0.6, -0.15))

    plt.savefig("Mutation_Network.png")

def one_mutation_away(seq1, seq2):
    mutations = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mutations+=1
        if mutations > 1:
            return False
    return mutations==1    


def get_intermediate_mutations(original, sequences):
    intermediate = []
    for s in sequences:
        for i, n in enumerate(s):
            if original[i] != n:
                inter = s[:i]+original[i]+s[i+1:]
                if inter not in intermediate:
                    intermediate.append(inter)
    return intermediate


bind_BA1 = BindingCalculator(known_to_neutralize="BA.1", weight_by_log_IC50=False)
bind_BA2 = BindingCalculator(known_to_neutralize="BA.2", weight_by_log_IC50=False)
bind_BA275 = BindingCalculator(known_to_neutralize="BA.2.75", weight_by_log_IC50=False)
bind_BA5 = BindingCalculator(known_to_neutralize="BA.5", weight_by_log_IC50=False)
bind_BQ11 = BindingCalculator(known_to_neutralize="BQ.1.1", weight_by_log_IC50=False)
bind_D614G = BindingCalculator(known_to_neutralize="D614G", weight_by_log_IC50=False)
bind_XBB = BindingCalculator(known_to_neutralize="XBB", weight_by_log_IC50=False)

def get_titer_table(sites):
    # testing against antibodies known to neutralize BA.1 BA.2 BA.275 BQ.1.1 D614G XBB
    columns = ["", "BA.1", "BA.2", "BA.2.75", "BA.5", "BQ.1.1", "D614G", "XBB"]
    titer = []
    for i in range(len(sites)):
        titer.append([
            "Seq"+str(i),
            bind_BA1.binding_retained({sites[i]}),
            bind_BA2.binding_retained({sites[i]}),
            bind_BA275.binding_retained({sites[i]}),
            bind_BA5.binding_retained({sites[i]}),
            bind_BQ11.binding_retained({sites[i]}),
            bind_D614G.binding_retained({sites[i]}),
            bind_XBB.binding_retained({sites[i]})
            ])
    frame = pd.DataFrame(columns=columns, data=titer)
    frame.to_csv("results/titer_table.csv", index=False)
    #titer.to_csv("results/titer_table.csv", index=False)


if __name__ == "__main__":
    data = pd.read_csv("results/fitness_change_of_antigenically_non_neutral_2_hop.csv", index_col=False)
    # Get the top value in the old sequence column
    old_sequence = data["old sequence"].iloc[0]
    #Down select
    data = data[["site", "fitness change", "antigen difference", "new sequence"]]
    #Select where anitgen difference above average
    mean = data["antigen difference"].mean()
    #Select only positive fitness change
    data = data[(data["fitness change"] > 0)]
    #Code for above average antigen difference, turned off right now
    #data = data[(data["antigen difference"] > mean)]
    data.to_csv("results/down_selected.csv", index=False)
    #List of sequences
    sequences = list(data["new sequence"])

    # get titer table values
    sites = data['site']
    get_titer_table(sites)

    original_length = len(sequences)

    # add original sequence to new sequences
    sequences.append(old_sequence)
    fitness_change = list(data["fitness change"])
    fitness_change.append(0)
    antigen_difference = list(data["antigen difference"])
    antigen_difference.append(0)

    # add intermediate mutations between
    # original sequences and new sequences
    intermediate_mutations = get_intermediate_mutations(old_sequence, sequences)
    sequences = sequences+intermediate_mutations
    fitness_change = fitness_change + ([0] * len(intermediate_mutations))
    antigen_difference = antigen_difference + ([0] * len(intermediate_mutations))

    #create_graph(sequences, fitness_change, antigen_difference, original_length)

    print("krewl")