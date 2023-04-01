import pandas as pd

if __name__ == "__main__":
    data = pd.read_csv("results/fitness_change_of_antigenically_non_neutral_2_hop.csv", index_col=False)
    # Get the top value in the old sequence column
    old_sequence = data["old sequence"].iloc[0]
    #Down select
    data = data[["fitness change", "antigen difference", "new sequence"]]
    #Select where anitgen difference above average
    mean = data["antigen difference"].mean()
    #Select only positive fitness change
    data = data[(data["fitness change"] > 0)]
    #Code for above average antigen difference, turned off right now
    #data = data[(data["antigen difference"] > mean)]
    data.to_csv("results/down_selected.csv", index=False)
    #List of sequences
    sequences = list(data["new sequence"])
    print("krewl")