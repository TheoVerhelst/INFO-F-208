from math import log
from copy import deepcopy

gap = "-"
aminoAcids = "ARNDCEQGHILKMFPSTWYV"

# The following probabilities has been gathered from the Swiss-prot database
# statistics page
p={ "E" : 0.0674, "H" : 0.0227, "S" : 0.066,  "K" : 0.0582, "Y" : 0.0292,
	"T" : 0.0535, "W" : 0.0109, "L" : 0.0965, "R" : 0.0553, "M" : 0.0241,
	"Q" : 0.0393, "A" : 0.0826, "V" : 0.0687, "F" : 0.0386, "C" : 0.0137,
	"N" : 0.0406, "P" : 0.0472, "I" : 0.0593, "D" : 0.0546, "G" : 0.0708}

logBase = 10

class Pssm:
    def __init__(self, filename):
        self.load(filename)
    
    def load(self, filename):
        sequences = []
        
        with open(filename) as file:
            for line in file:
                line = line.strip()
                if not line.startswith(">"):
                    sequences.append(line.upper())
        
        self.numberColumns = len(sequences[0])
        matrix = [{aa : 0 for aa in aminoAcids} for _ in range(self.numberColumns)]
        beta = len(sequences) ** 0.5
        self.matrix = deepcopy(matrix)
        
        for sequence in sequences:
            for i, aa in enumerate(sequence):
                if aa in aminoAcids:
                    matrix[i][aa] += 1
                    
        
        for i, column in enumerate(matrix):
            alpha = sum(column.values())
            for aa, count in column.items():
                frequency = count / len(sequences)
                probability = (alpha * frequency + beta * p[aa]) / (alpha + beta)
                self.matrix[i][aa] = log(probability / p[aa], logBase)
    
    def __getitem__(self, column, aa):
        return self.matrix[column][aa]
    
    def __str__(self):
        res = "  | "
        cellWidth = 5
        for i in range(self.numberColumns):
            res += str(i).rjust(cellWidth) + " | "
        res += "\n  " + "+-------" * self.numberColumns + "\n"
        for aa in aminoAcids:
            res += aa + " | "
            for column in self.matrix:
                res += str(column[aa])[:cellWidth].rjust(cellWidth) + " | "
            res += "\n"
        res += "\n"
        return res

def main():
    muscle = "msaresults-MUSCLE.fasta"
    clustal = "msaresults-CLUSTAL.fasta"
    #pssm = Pssm(muscle)
    pssm = Pssm(clustal)
    print(pssm)
    
if __name__ == "__main__":
    main()
