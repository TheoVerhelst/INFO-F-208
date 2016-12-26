from math import log, sqrt
import pickle
from matplotlib import pyplot
from copy import deepcopy

amino_acids = "ARNDCEQGHILKMFPSTWYVBZXJOU"
structures = "HCTE"
window_width = 8
log_base = 2
terminal_width = 140

def parse_cath(filename):
    res = []
    with open(filename) as file:
        for line in file:
            res.append({"file" : line[0:4], "sequence" : line[4]})
    return res

def parse_dssp(filename, sequence):
    aa_sequence = ""
    struct_sequence = ""
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
            
            if identifier.lower() == sequence.lower():
                amino_acid = line[amino_acid_column].upper()
                structure = line[structure_column].upper()
                
                if structure not in classes.keys():
                    raise ValueError("File " + filename + ", line " + str(i) + ": Unknown secondary structure \"" + structure + "\"")
                if amino_acid not in amino_acids:
                    raise ValueError("File " + filename + ", line " + str(i) + ": Unknown amino acid \"" + amino_acid + "\"")
            
                aa_sequence += amino_acid
                struct_sequence += classes[structure]
                
        return (aa_sequence, struct_sequence)
        
def get_all_sequences(dssp_directory, cath_info_filename):
    res = []
    for parse_data in parse_cath(cath_info_filename):
        res.append(parse_dssp(dssp_directory + parse_data["file"] + ".dssp", parse_data["sequence"]))
    return res

def train(dssp_directory, cath_info_filename):
    sequences = get_all_sequences(dssp_directory, cath_info_filename)
    f_s = {s : 0 for s in structures}
    f_s_r = {(aa, s) : 0 for s in structures for aa in amino_acids}
    f_s_r_rm = {(aa, s, aa_i) : 0
                            for aa in amino_acids
                            for s in structures
                            for aa_i in amino_acids}
                            
    print("Training on", len(sequences), "sequences")
                            
    for aa_sequence, struct_sequence in sequences:
        for i, (aa, struct) in enumerate(zip(aa_sequence, struct_sequence)):
            f_s[struct] += 1
            f_s_r[(aa, struct)] += 1
            for aa_m, struct_m in sequence[max(0, i - window_width) : i] \
                                + sequence[i + 1 : i + window_width + 1]:
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

def print_prediction(aa_sequence, real_sequence, predicted_sequence):
    predicted_sequence = "".join(predicted["max"] for predicted in predicted_sequence)
    first_line = "Amino acids sequence     :"
    
    index = terminal_width - len(first_line) - 1
    print(first_line, aa_sequence[:index])
    print("Real structure sequence  :", real_sequence[:index])
    print("Predicted struct sequence:", predicted_sequence[:index])
    print()
    
    while index < len(aa_sequence):
        print(aa_sequence[index:index + terminal_width])
        print(real_sequence[index:index + terminal_width])
        print(predicted_sequence[index:index + terminal_width])
        print()
        index += terminal_width

def run_tests(dssp_directory, cath_info_filename, f_s, f_s_r, f_s_r_rm, plot=False):
    sequences = get_all_sequences(dssp_directory, cath_info_filename)
    results = []
    
    for i, (aa_sequence, real_sequence) in enumerate(sequences):
        predicted_sequence = predict(f_s, f_s_r, f_s_r_rm, aa_sequence)
        print("Sequence", i + 1)
        print_prediction(aa_sequence, real_sequence, predicted_sequence)
        zipped = list(zip(real_sequence, predicted_sequence))
        
        Q3 =  sum(1 for pos in zipped if pos[0] == pos[1]["max"]) / len(zipped)
        MCC = {}
        
        for structure in structures:
            TP, TN, FP, FN = 0, 0, 0, 0
            P = sum(1 for pos in zipped if structure == pos[1]["max"])
            N = len(zipped) - P
            roc_curve = []
            # Sort by score from this structure, in order to construct the ROC curve
            zipped.sort(key=lambda pos : pos[1][structure], reverse=True)
            
            for real, predicted in zipped:
                roc_curve.append(((TN + FN)/N, (TP + FP)/P))
                    
                if predicted["max"] == structure:
                    if real == predicted["max"]:
                        TP += 1
                    else:
                        FP += 1
                else:
                    if real != structure:
                        TN += 1
                    else:
                        FN += 1

            try:
                MCC[structure] =  (TP * TN - FP * FN) / \
                        sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
            except ZeroDivisionError:
                MCC[structure] = float("nan")
            
            if plot:
                roc_curve.append(((TN + FN)/N, (TP + FP)/P))
                pyplot.plot(*zip(*roc_curve), label="Structure " + structure)
                pyplot.ylabel("FPR")
                pyplot.ylabel("TPR")
        
        print("Q3 =", Q3, ", MCC =", MCC)
        pyplot.title("ROC curve")
        pyplot.legend()
        pyplot.show()
        print()
        print("*" * terminal_width)
        print()
        

def main():
    # Script parameters
    already_pickled = True
    pickle_filename = "frequencies.pickle"
    dssp_directory = "dataset/dssp/"
    cath_info_filename = "dataset/CATH_info.txt"
    dssp_directory_test = "dataset/dssp_test/"
    cath_info_filename_test = "dataset/CATH_info_test.txt"
    
    print("Gathering statistics about the training set...")
    if already_pickled:
        with open(pickle_filename, 'rb') as file:
            frequencies = pickle.load(file)
    else:
        frequencies = train(dssp_directory, cath_info_filename)
        with open(pickle_filename, 'wb') as file:
            pickle.dump(frequencies, file)
            
    print("Running tests...")
    run_tests(dssp_directory_test, cath_info_filename_test, *frequencies, plot=True)

if __name__ == "__main__":
    main()
