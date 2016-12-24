from math import log, sqrt
import pickle
from matplotlib import pyplot
from copy import deepcopy

amino_acids = "ARNDCEQGHILKMFPSTWYVBZXJOU"
structures = "HCTE"
window_width = 8
log_base = 10

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
        " " : "C"
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
                    raise ValueError("File " + filename + ", line " + str(i) + ": Unknown secondary structure \"" + structure + "\"")
                if amino_acid not in amino_acids:
                    raise ValueError("File " + filename + ", line " + str(i) + ": Unknown amino acid \"" + amino_acid + "\"")
                res.append((amino_acid, classes[structure]))
        return res
        
def get_all_sequences(dssp_directory, cath_info_filename):
    res = []
    for parse_data in parse_cath(cath_info_filename):
        res.append(parse_dssp(dssp_directory + parse_data["file"] + ".dssp", parse_data["sequence"]))
    return res

def train(dssp_directory, cath_info_filename):
    sequences = get_all_sequences(dssp_directory, cath_info_filename)
    print("Training on", len(sequences), "sequences")
    
    f_s = {s : 0 for s in structures}
    f_s_r = {(aa, s) : 0 for s in structures for aa in amino_acids}
    f_s_r_rm = {(aa, s, aa_i) : 0
                            for aa in amino_acids
                            for s in structures
                            for aa_i in amino_acids}
    for i, sequence in enumerate(sequences):
        for i, (aa, struct) in enumerate(sequence):
            f_s[struct] += 1
            f_s_r[(aa, struct)] += 1
            for aa_m, struct_m in sequence[max(0, i - window_width) : i] + sequence[i + 1 : i + window_width + 1]:
                f_s_r_rm[(aa, struct, aa_m)] += 1
    return f_s, f_s_r, f_s_r_rm

def predict(f_s, f_s_r, f_s_r_rm, sequence):
    res = []
    for i, aa in enumerate(sequence):
        scores = {}
        for s in structures:
            n_s = [s2 for s2 in structures if s2 != s]
            score = log(f_s_r[(aa, s)] / sum(f_s_r[(aa, s2)] for s2 in n_s), log_base)
            score += log(sum(f_s[s2] for s2 in n_s) / f_s[s], log_base)
            for aa_m in sequence[max(0, i - window_width) : i] + sequence[i + 1 : i + window_width + 1]:
                score += log(f_s_r_rm[(aa, s, aa_m)]
                        / sum(f_s_r_rm[(aa, s2, aa_m)] for s2 in n_s), log_base)
                score += log(sum(f_s_r[(aa, s2)] for s2 in n_s) / f_s_r[(aa, s)], log_base)
            scores[s] = score
        scores["max"] = max(scores, key=lambda k:scores[k])
        res.append(scores)
    return res
                    
def run_tests(dssp_directory, cath_info_filename, f_s, f_s_r, f_s_r_rm):
    sequences = get_all_sequences(dssp_directory, cath_info_filename)
    """
    big_sequence = []
    for sequence in sequences:
        big_sequence += sequence
    sequences = [big_sequence]
    """
    
    for i, sequence in enumerate(sequences):
        print("Predicting sequence " + str(i + 1) + "...")
        aa_sequence, real_sequence = zip(*sequence)
        predicted_sequence = predict(f_s, f_s_r, f_s_r_rm, aa_sequence)
        seq_data = [{"real" : real,
                    "predicted" : predicted}
                    # The zip expression unzips `sequence`, and zips the result with `predicted_sequence`
                    for real, predicted in zip(real_sequence, predicted_sequence)]
        
        Q3 =  sum(1 for pos_data in seq_data if pos_data["real"] == pos_data["predicted"]["max"]) \
                / len(seq_data)
        MCC = {}
        
        for structure in structures:
            TP, TN, FP, FN = 0, 0, 0, 0
            P = sum(1 for pos_data in seq_data if structure == pos_data["predicted"]["max"])
            N = len(seq_data) - P
            TP_roc = 0
            FP_roc = 0
            prev_score = float("-inf")
            roc_curve_x, roc_curve_y = [], []
            # Sort by score from this structure, in order to construct the ROC curve
            seq_data.sort(key=lambda pos_data : pos_data["predicted"][structure], reverse=True)
            
            for pos_data in seq_data:
                score = pos_data["predicted"][structure]
                
                if score != prev_score:
                    roc_curve_x.append(FP_roc/N)
                    roc_curve_y.append(TP_roc/P)
                    
                if pos_data["predicted"]["max"] == structure:
                    if pos_data["real"] == pos_data["predicted"]["max"]:
                        TP += 1
                    else:
                        FP += 1
                    TP_roc += 1
                else:
                    if pos_data["real"] != structure:
                        TN += 1
                    else:
                        FN += 1
                    FP_roc += 1

            try:
                MCC[structure] =  (TP * TN - FP * FN) / \
                        sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
            except ZeroDivisionError:
                MCC[structure] = float("nan")
            
            roc_curve_x.append(FP_roc/N)
            roc_curve_y.append(TP_roc/P)
            pyplot.plot(roc_curve_x, roc_curve_y)
            pyplot.ylabel("FPR")
            pyplot.ylabel("TPR")
            pyplot.title("ROC curve for " + structure)
            pyplot.show()
        
        print("Q3 =", Q3, ", MCC =", MCC)

def main():
    # Script parameters
    already_pickled = True
    pickle_filename = "frequencies.pickle"
    dssp_directory = "dataset/dssp/"
    cath_info_filename = "dataset/CATH_info.txt"
    dssp_directory_test = "dataset/dssp_test/"
    cath_info_filename_test = "dataset/CATH_info_test.txt"
    
    if already_pickled:
        with open(pickle_filename, 'rb') as file:
            frequencies = pickle.load(file)
    else:
        frequencies = train(dssp_directory, cath_info_filename)
        with open(pickle_filename, 'wb') as file:
            pickle.dump(frequencies, file)
            
    run_tests(dssp_directory_test, cath_info_filename_test, *frequencies)

if __name__ == "__main__":
    main()
