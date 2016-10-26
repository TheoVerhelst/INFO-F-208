from enum import Enum
from copy import deepcopy

AminoAcid = Enum("AminoAcid", "A R N D C E Q G H I L K M F P S T W Y V B Z X")
# Make the str() function print only the amino acid letter by redefining __str__
AminoAcid.__str__ = lambda self: self.name[-1]

class Sequence:
    """Represents a sequence of amino acids (with no further semantic).
    It is suited to compute the best alignment between two sequences."""
    
    def __init__(self, filename, sequenceId):
        """Constructor. Loads the sequence from a file.
        
        Parameters:
            filename: the filename of the file containing the sequence.
            sequenceId: the identifier of the sequence in the file.
                This is the string between the two first bars "|" on the line
                preceding the sequence.
            """
        self.loadFromFile(filename, sequenceId)
    
    def __getitem__(self, index):
        """Returns the elements at the specified index. Allows the syntax
        sequence[i].
        """
        return self.sequence[index]
        
    def __len__(self):
        """Returns the lenght of the sequence. Allows the syntax len(sequence)."""
        return len(self.sequence)

    def loadFromFile(self, filename, sequenceId):
        """Loads the sequence data from a file.
        
        Parameters:
            filename: the filename of the file containing the sequence.
            sequenceId: the identifier of the sequence in the file.
                This is the string between the two first bars "|" on the line
                preceding the sequence.
        """
        foundSequence = False
        sequence = ""
        with open(filename) as file:
            for line in file:
                if not foundSequence:
                    # We found the sequence if this is the right identifier line
                    foundSequence = (line[0] == ">" and line.split("|")[1] == sequenceId)
                else:
                    if line[0] != ">":
                        sequence += line
                    else:
                        # Stop adding lines to the sequence
                        foundSequence = False
                        
        if sequence == "":
            raise RuntimeError("Sequence not found in the specified file.")

        try:
            sequence = [AminoAcid[character] for character in sequence if not character.isspace()]
        except KeyError as e:
            raise RuntimeError("An amino acid is not valid in the sequence: " + str(e))

        # We didn't worked on self.sequence directly to preserve the internal state in case of failure
        self.sequence = sequence
    
class Score:
    """Class used to compute the score of substitution between two amino acids,
    by using well known matrix such as BLOSUM.
    """
    
    def __init__(self, filename):
        """Contructor. Loads the matrix from a file.
        
        Parameter:
            filename: the filename of the matrix to load. It must be a .iij file
            as given in the blosum archive.
        """
        self.loadFromFile(filename)
    
    def getScore(self, aminoAcidA, aminoAcidB):
        """Returns the score of matching between two amino acids according to
        the loaded matrix.
        """
        indexA = self.indices.index(aminoAcidA)
        indexB = self.indices.index(aminoAcidB)
        
        # The matrix is triangular, so we need to test the two possibilities for indexing
        # since one of the two would result in an index error
        if indexA > indexB:
            return self.matrix[indexA][indexB]
        else:
            return self.matrix[indexB][indexA]
    
    def loadFromFile(self, filename):
        """Loads a substitution matrix from a file.
        
        Parameter:
            filename: the filename of the matrix to load. It must be a .iij file
            as given in the blosum archive.
        """
        self.matrix = []
        self.indices = []
        
        with open(filename) as file:
            foundAminoAcidList = False
            
            for line in file:
                # Skip commentary lines
                if line.strip().startswith("#"):
                    continue
                
                if not foundAminoAcidList:
                    # The first non-commentary line is the list of amino acids
                    # We reload this list because some files may order the acids differently
                    self.indices = [AminoAcid[character] for character in line.split()]
                    foundAminoAcidList = True
                else:
                    self.matrix.append([int(number) for number in line.split()])

def needlmanWunsch(score, sequence1, sequence2, openPenalty, extendPenalty):
    """Applies de Needleman-Wunsch algorithm, given a scoring matrix,
    two sequences and a gap weight.
    
    Parameters:
        score: an instance of Score
        sequence1, sequence2: two instances of Sequence
        g: the weight of a gap
    
    Return value: The set of solutions for aligning these two sequences.
    """
    m, n = len(sequence1) + 1, len(sequence2) + 1
    
    # Construct the alignment matrix first
    S = [[0] * n for i in range(m)]
    backtrace = [["0"] * n for i in range(m)]
    
    V, W = deepcopy(S), deepcopy(S)
        
    for i in range(1, m):
        for j in range(1, n):
            V[i][j] = max(S[i - 1][j] + openPenalty, V[i - 1][j] + extendPenalty, 0)
            W[i][j] = max(S[i][j - 1] + openPenalty, W[i][j - 1] + extendPenalty, 0)
            
            # We use a dict to find the value and the origin string for the
            # current cell at the same time
            choices = {
                "0" : 0, # For the local alignment
                "T" : V[i][j],
                "L" : W[i][j],
                # We decrease the index for accessing the sequences because the
                # sequences have one less elements than the matrices
                "D" : S[i - 1][j - 1] + score.getScore(sequence1[i - 1], sequence2[j - 1])
            }
            # Find the maximum value in the dict
            S[i][j] = max(choices.values())
            # Add all origins that can lead to this maximum value in the backtrace
            backtrace[i][j] = "".join([key for key in choices if choices[key] == S[i][j]])

    # Then find the alignment from the backtrace, by starting from the maximum in S
    return findAlignment(backtrace, *findMaximum(S))

def findMaximum(matrix):
    """Finds one maximum in a matrix.
    
    Parameters:
        matrix: the matrix
    
    Return value: a tuple containing the (i, j) coordinates of the maximum
    """
    maximum = float("-inf")
    for i, line in enumerate(matrix):
        for j, value in enumerate(line):
            if value > maximum:
                maximum = value
                imax, jmax = i, j
    return imax, jmax
    
# TODO: remove this function
def printMatrix(matrix, sequence1, sequence2):
    size = 8
    print(" ", "".join([" " * size] + [str(s).rjust(size) for s in sequence2]))
    for i, line in enumerate(matrix):
        if i > 0:
            print(sequence1[i - 1], end=" ")
        else:
            print("  ", end="")
        for cell in line:
            print(str(cell).rjust(size), end='')
        print()
    
def findAlignment(backtrace, i, j, currentSolution = "", solutions = None):
    """Recursively iterates on the backtrace matrix to find the path from the
    bottom-right cell to the top-left cell.
    
    Parameters:
        backtrace: the matrix indicating the origin of each cell
        i, j: the coordinates of the currently explored cell
        currentSolution: a string reprensenting the currently explored solution
        solution: the list of all valid solution found so far
    
    Return value: the list of solutions
    """
    if solutions == None:
        solutions = []
    
    # Stop if we are in the top-left cell, or in a zero cell
    if (i == 0 and j == 0) or "0" in backtrace[i][j]:
        solutions.append({"origin" : (i, j), "path" : currentSolution})
    else:
        for possibility in backtrace[i][j]:
            currentSolution += possibility
            if possibility == "T":
                findAlignment(backtrace, i - 1, j, currentSolution, solutions)
            elif possibility == "L":
                findAlignment(backtrace, i, j - 1, currentSolution, solutions)
            elif possibility == "D":
                findAlignment(backtrace, i - 1, j - 1, currentSolution, solutions)
            currentSolution = currentSolution[:-1]
    return solutions
    
def printAlignment(sequence1, sequence2, solutions):
    """Prints to stdout all the alignemnts given in parameters.
    
    Parameters:
        sequence1, sequence2: the two aligned instances of Sequence
        solutions: the list of all possible alignments
    """
    gapChar = "-"
    indelChar = " "
    conservationChar = ":"
    mutationChar = "."
    maxLineLength = 120
    print(len(solutions), "solutions were found.")
    
    for k, solution in enumerate(solutions):
        sequence1Str, sequence2Str, midStr= "", "", ""
        i, j = solution["origin"]
        i -= 1
        j -= 1
        
        # Iterate over the solution backward
        for origin in reversed(solution["path"]):
            if origin == "T":
                i += 1
                sequence1Str += str(sequence1[i])
                sequence2Str += gapChar
                midStr += indelChar
            elif origin == "L":
                j += 1
                sequence1Str += gapChar
                sequence2Str += str(sequence2[j])
                midStr += indelChar
            elif origin == "D":
                i += 1
                j += 1
                sequence1Str += str(sequence1[i])
                sequence2Str += str(sequence2[j])
                midStr += conservationChar if sequence1[i] == sequence2[j] else mutationChar
                
        print("Solution No " + str(k + 1) + ":")
        for k in range(0, ((len(sequence1Str) - 1) // maxLineLength) + 1):
            slice = (k * maxLineLength, min((k + 1) * maxLineLength, len(sequence1Str)))
            print("Sequence 1:", sequence1Str[slice[0]:slice[1]])
            print("           ", midStr[slice[0]:slice[1]])
            print("Sequence 2:", sequence2Str[slice[0]:slice[1]])
            print()
        print()
        
if __name__ == "__main__":
    # Script parameters
    scoringMatrixFilename = "blosum/blosum50.iij"
    sequence1File = "maguk-sequences.fasta"
    sequence1Id = "Q12959"
    sequence2File = "maguk-sequences.fasta"
    sequence2Id = "Q92796"
    openPenalty = -12
    extendPenalty = -2

    score = Score(scoringMatrixFilename)
    sequence1 = Sequence(sequence1File, sequence1Id)
    sequence2 = Sequence(sequence2File, sequence2Id)

    printAlignment(sequence1, sequence2, needlmanWunsch(score, sequence1, sequence2, openPenalty, extendPenalty))
