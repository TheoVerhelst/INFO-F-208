from enum import Enum

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
                This is the string between the two first bars "|" on the line preceding the sequence.
            """
        self.loadFromFile(filename, sequenceId)
    
    def __getitem__(self, index):
        """Returns the elements at the specified index. Allows the syntax sequence[i]."""
        return self.sequence[index]
        
    def __len__(self):
        """Returns the lenght of the sequence. Allows the syntax len(sequence)."""
        return len(self.sequence)

    def loadFromFile(self, filename, sequenceId):
        """Loads the sequence data from a file.
        
        Parameters:
            filename: the filename of the file containing the sequence.
            sequenceId: the identifier of the sequence in the file.
                This is the string between the two first bars "|" on the line preceding the sequence.
        """
        found = False
        with open(filename) as file:
            for line in file:
                # If we found the needed identifier line
                if line[0] == ">" and line.split("|")[1] == sequenceId:
                    found = True
                    sequence = file.readline()
                    break

        if not found:
            raise RuntimeError("Sequence not found in the specified file.")

        try:
            sequence = [AminoAcid[character] for character in sequence if not character.isspace()]
        except KeyError as e:
            raise RuntimeError("An amino acid is not valid in the sequence: " + str(e))

        # We didn't worked on self.sequence directly to preserve the internal state in case of failure
        self.sequence = sequence
    
class Score:
    """Class used to compute the score of substitution between two amino acids, by using well known matrix
    such as the BLOSUM."""
    
    def __init__(self, filename):
        """Contructor. Loads the matrix from a file.
        
        Parameter:
            filename: the filename of the matrix to load. It must be a .iij file as given in the blosum archive.
        """
        self.loadFromFile(filename)
    
    def getScore(self, aminoAcidA, aminoAcidB):
        """Returns the score of matching between two amino acids according to the loaded matrix."""
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
            filename: the filename of the matrix to load. It must be a .iij file as given in the blosum archive.
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
                    self.matrix.append([float(number) for number in line.split()])

def needlmanWunsch(score, sequence1, sequence2, g):
    m, n = len(sequence1) + 1, len(sequence2) + 1
    
    # Construct the alignment matrix first
    matrix = [[0] * n for i in range(m)]
    backtrace = [["0"] * n for i in range(m)]
    
    for i in range(1, m):
        matrix[i][0] = i * g
        backtrace[i][0] = "T"
        
    for j in range(1, n):
        matrix[0][j] = j * g
        backtrace[0][j] = "L"
        
    for i in range(1, m):
        for j in range(1, n):
            # We use a dict to find the value and the origin string for the
            # current cell at the same time
            choices = {
                "T" : matrix[i - 1][j] + g,
                "L" : matrix[i][j - 1] + g,
                "D" : matrix[i - 1][j - 1] + score.getScore(sequence1[i - 1], sequence2[j - 1])
            }
            # Find the maximum value in the dict
            maximumValue = max(choices.values())
            matrix[i][j] = maximumValue
            # Add all origins that can lead to this maximum value in the backtrace
            backtrace[i][j] = "".join([key for key in choices if choices[key] == maximumValue])

    # Then find the alignment from the backtrace
    return findAlignment(backtrace, m - 1, n - 1)
    
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
    if solutions == None:
        solutions = []
    
    if i == 0 and j == 0:
        solutions.append(currentSolution)
    else:
        for possibility in backtrace[i][j]:
            if possibility == "T":
                findAlignment(backtrace, i - 1, j, currentSolution + possibility, solutions)
            elif possibility == "L":
                findAlignment(backtrace, i, j - 1, currentSolution + possibility, solutions)
            elif possibility == "D":
                findAlignment(backtrace, i - 1, j - 1, currentSolution + possibility, solutions)
    return solutions
    
def printAlignment(sequence1, sequence2, solutions):
    gapChar = "-"
    indelChar = " "
    conservationChar = ":"
    mutationChar = "."
    print(len(solutions), "solutions were found.")
    for k, solution in enumerate(solutions):
        sequence1Str, sequence2Str, midStr= "", "", ""
        i, j = -1, -1
        
        # Iterate over the solution backward
        for origin in reversed(solution):
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
        print("Sequence 1: ", sequence1Str)
        print("            ", midStr)
        print("Sequence 2: ", sequence2Str)
        print()
        
# Script parameters
scoringMatrixFilename = "blosum/blosum80.iij"
sequence1File = "SH3-sequence.fasta"
sequence1Id = "P12931"
sequence2File = "SH3-sequence.fasta"
sequence2Id = "P62993"
g = -4

score = Score(scoringMatrixFilename)
sequence1 = Sequence(sequence1File, sequence1Id)
sequence2 = Sequence(sequence2File, sequence2Id)

printAlignment(sequence1, sequence2, needlmanWunsch(score, sequence1, sequence2, g))
