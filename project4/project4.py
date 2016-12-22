from math import log, sqrt
import pickle
import matplotlib.pyplot as plt
from copy import deepcopy

amino_acids = "ARNDCEQGHILKMFPSTWYVBZXJUO"
structures = "HCTE"
window_width = 8
log_base = 2

def parse_cath(filename):
    res = []
    with open(filename) as file:
        for line in file:
            res.append({"file" : line[0:4], "sequence" : line[4]})
    return res

def parse_dssp(filename, sequence):
    res = []
    classes = \
    {
        "H" : "H",
        "G" : "H",
        "I" : "H",
        "E" : "E",
        "B" : "E",
        "T" : "T",
        "C" : "C",
        "S" : "C",
        " "  : "C"
    }
    # These values has been found manually, it seems that all the database
    # matches this format
    identifier_column = 11
    amino_acid_column = 13
    structure_column = 16
    
    with open(filename) as file:
        line = ""
        i = 0
        
        while not line.startswith("#"):
            i += 1
            line = file.readline().strip()
            
        for line in file:
            i += 1
            identifier = line[identifier_column].upper()
            amino_acid = line[amino_acid_column].upper()
            structure = line[structure_column].upper()
            if identifier.lower() == sequence.lower():
                # If there is no data for secondary structure, then the split
                # returned the next column, so we force to have empty string
                if structure not in classes.keys():
                    raise Value_error("File " + filename + ", line " + str(i) + ": Unknown secondary structure \"" + structure + "\"")
                if amino_acid not in amino_acids:
                    raise Value_error("File " + filename + ", line " + str(i) + ": Unknown amino acid \"" + amino_acid + "\"")
                res.append((amino_acid, classes[structure]))
        return res
        
def get_all_sequences(dssp_prefix, cath_info_filename):
    res = []
    for parse_data in parse_cath(cath_info_filename):
        res.append(parse_dssp(dssp_prefix + parse_data["file"] + ".dssp", parse_data["sequence"]))
    return res

def train():
    dssp_prefix = "dataset/dssp/"
    cath_info_filename = "dataset/CATH_info.txt"
    sequences = get_all_sequences(dssp_prefix, cath_info_filename)
    
    f_s = {s : 0 for s in structures}
    f_s_r = {(aa, s) : 0 for s in structures for aa in amino_acids}
    f_s_r_rm = {(aa, s, aa_i) : 0
                            for aa in amino_acids
                            for s in structures
                            for aa_i in amino_acids
                            for i in range(-window_width, window_width + 1)
                            if i != 0}
    for sequence in sequences:
        for i, (aa, struct) in enumerate(sequence):
            f_s[struct] += 1
            f_s_r[(aa, struct)] += 1
            for m in range(-window_width, window_width + 1):
                if m != 0 and i + m >= 0 and i + m < len(sequence):
                    f_s_r_rm[(aa, struct, sequence[i + m][0])] += 1
    
    return f_s, f_s_r, f_s_r_rm

def predict(f_s, f_s_r, f_s_r_rm, sequence):
    res = []
    for i, aa in enumerate(sequence):
        scores = {}
        for s in structures:
            other_s = [sp for sp in structures if sp != s]
            score = log(f_s_r[(aa, s)] / sum(f_s_r[aa, sp] for sp in other_s), log_base)
            score += log(sum(f_s[sp] for sp in other_s) / f_s[s])
            for m in range(-window_width, window_width + 1):
                if m != 0 and i + m >= 0 and i + m < len(sequence):
                    score += log(f_s_r_rm[(aa, s, sequence[i+m])]
                            / sum(f_s_r_rm[(aa, sp, sequence[i+m])] for sp in other_s), log_base)
                    score += log(sum(f_s_r[aa, sp] for sp in other_s) / f_s_r[(aa, s)], log_base)
            scores[s] = score
        #scores["max"] = max(scores, key=lambda k:scores[k])
        res.append({"max" : max(scores, key=lambda k:scores[k])})
    return res
                    
def run_tests(f_s, f_s_r, f_s_r_rm):
    dssp_prefix = "dataset/dssp_test/"
    cath_info_filename = "dataset/CATH_info_test.txt"
    sequences = get_all_sequences(dssp_prefix, cath_info_filename)
    
    for i, sequence in enumerate(sequences):
        print("Sequence", i + 1)
        real_sequence = [t[1] for t in sequence]
        predicted_sequence = predict(f_s, f_s_r, f_s_r_rm, real_sequence)
        seq_data = [{"real" : real,
                    "predicted" : predicted}
                    # The zip expression unzips `sequence`, and zips the result with `predicted_sequence`
                    for real, predicted in zip(real_sequence, predicted_sequence)]
        for pos_data in seq_data:
            if pos_data["real"] == pos_data["predicted"]["max"]:
                print("X", end="")
            else:
                print(" ", end="")
        print()
        
        Q3 =  sum(1 for pos_data in seq_data if pos_data["real"] == pos_data["predicted"]["max"]) \
                / len(seq_data)
        MCC = {}
        
        for structure in structures:
            TP, TN, FP, FN = 0, 0, 0, 0
            # Sort by score from this structure, in order to construct the ROC curve
            #seq_data.sort(key=lambda pos_data : pos_data["predicted"][structure])
            for pos_data in seq_data:
                if pos_data["predicted"]["max"] == structure:
                    if pos_data["real"] == pos_data["predicted"]["max"]:
                        TP += 1
                    else:
                        FP += 1
                else:
                    if pos_data["real"] == pos_data["predicted"]["max"]:
                        TN += 1
                    else:
                        FN += 1
            try:
                MCC[structure] =  (TP * TN - FP * FN) / \
                        sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
            except ZeroDivisionError:
                MCC[structure] = float("nan")
    
        print("Q3 =", Q3, ", MCC =", MCC)

def main():
    already_pickled = True
    pickle_filename = "frequencies.pickle"
    if already_pickled:
        with open(pickle_filename, 'rb') as file:
            frequencies = pickle.load(file)
    else:
        frequencies = train()
        with open(pickle_filename, 'wb') as file:
            pickle.dump(frequencies, file)
            
    run_tests(*frequencies)

if __name__ == "__main__":
    main()
