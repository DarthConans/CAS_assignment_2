import pandas as pd

def to_number(val):
    val = val[1:-1]
    ret = None
    while ret is None:
        try:
            ret = int(val)
        except:
            val = val[:-1]
    return ret



threshold = 0


fit = pd.read_csv("data/aamut_fitness_all.csv")
fit["number"] = fit["mutation"].apply(to_number)
nans = fit[(fit["aa_fitness"] != fit["aa_fitness"]) | (fit["aa_fitness"] < threshold)]
non_nas = fit[(fit["aa_fitness"] == fit["aa_fitness"]) & (fit["aa_fitness"] > threshold)]
increase_fitness = non_nas.shape[0] / fit.shape[0] * 100.0
do_not_increase = nans.shape[0] / fit.shape[0] * 100.0

print(f"{increase_fitness}% INCREASE FITNESS")
print(f"{do_not_increase}% DO NOT INCREASE FITNESS")
print("krewl")