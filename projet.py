from enum import Enum
from copy import deepcopy
import cProfile

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

class Aligner:
    def __init__(self, score, sequence1, sequence2, openPenalty, extendPenalty):
        self.score = score
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.openPenalty = openPenalty
        self.extendPenalty = extendPenalty
        self.solutions = []
        # Execute the algorithm
        self.needlmanWunsch()
        self.printAlignment()
    
    def needlmanWunsch(self):
        """Applies de Needleman-Wunsch algorithm, given a scoring matrix,
        two sequences and a gap weight.
        
        Parameters:
            score: an instance of Score
            sequence1, sequence2: two instances of Sequence
            g: the weight of a gap
        
        Return value: The set of solutions for aligning these two sequences.
        """
        m, n = len(self.sequence1) + 1, len(self.sequence2) + 1
        
        # Construct the alignment matrix first
        self.S = [[0] * n for i in range(m)]
        self.backtrace = [["0"] * n for i in range(m)]
        
        V, W = deepcopy(self.S), deepcopy(self.S)
            
        for i in range(1, m):
            for j in range(1, n):
                V[i][j] = max(self.S[i - 1][j] + self.openPenalty, V[i - 1][j] + self.extendPenalty, 0)
                W[i][j] = max(self.S[i][j - 1] + self.openPenalty, W[i][j - 1] + self.extendPenalty, 0)
                
                # We use a dict to find the value and the origin string for the
                # current cell at the same time
                choices = {
                    "0" : 0, # For the local alignment
                    "T" : V[i][j],
                    "L" : W[i][j],
                    # We decrease the index for accessing the sequences because the
                    # sequences have one less elements than the matrices
                    "D" : self.S[i - 1][j - 1] + self.score.getScore(self.sequence1[i - 1], self.sequence2[j - 1])
                }
                # Find the maximum value in the dict
                self.S[i][j] = max(choices.values())
                # Add all origins that can lead to this maximum value in the backtrace
                self.backtrace[i][j] = "".join([key for key in choices if choices[key] == self.S[i][j]])

        # Find the alignment, starting from the cell with maximum value in S
        self.findAlignment(*self.findMaximum(self.S))

    def findMaximum(self, matrix):
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
        
    #TODO remove this
    def printMatrix(self, matrix):
        """Prints a matrix.
        This function is used to debug an test, it is not part of the final
        output.
        """
        cellSize = 8
        print(" " * (cellSize + 1), "".join([str(aa).rjust(cellSize) for aa in self.sequence2]))
        for i, line in enumerate(matrix):
            print(self.sequence1[i - 1] if i > 0 else " ", end=" ")
            for cell in line:
                print(str(cell).rjust(cellSize), end="")
            print()
        
    def findAlignment(self, i, j, currentSolution = ""):
        """Recursively iterates on the backtrace matrix to find valid paths from
        the cell at position (i, j) to a cell with value 0.
        All paths meeting this criteria are appended to self.solutions.
        
        Parameters:
            i, j: the coordinates of the currently explored cell
            currentSolution: a string reprensenting the currently explored solution
        """
        
        # Stop if we are in the top-left cell, or in a zero cell
        if (i == 0 and j == 0) or "0" in self.backtrace[i][j]:
            self.solutions.append({"origin" : (i, j), "path" : currentSolution})
        else:
            for possibility in self.backtrace[i][j]:
                currentSolution += possibility
                if possibility == "T":
                    self.findAlignment(i - 1, j, currentSolution)
                elif possibility == "L":
                    self.findAlignment(i, j - 1, currentSolution)
                elif possibility == "D":
                    self.findAlignment(i - 1, j - 1, currentSolution)
                currentSolution = currentSolution[:-1]
        
    def printAlignment(self):
        """Prints all the alignments found so far in self.solutions.
        """
        gapChar = "-"
        indelChar = " "
        conservationChar = ":"
        mutationChar = "."
        # The maximum number of amino acid displayed on a line
        maxLineLength = 120
        # The interval between two index hints dipslayed above or below the sequence
        indexHintInterval = 10
        
        print(len(self.solutions), "solutions were found.")
        
        for k, solution in enumerate(self.solutions):
            sequence1Str, sequence2Str, midStr= "", "", ""
            indexHints1, indexHints2 = "", ""
            i, j = solution["origin"]
            i -= 1
            j -= 1
            
            # Iterate over the solution backward
            for origin in reversed(solution["path"]):
                if origin == "T":
                    i += 1
                    sequence1Str += str(self.sequence1[i])
                    sequence2Str += gapChar
                    midStr += indelChar
                elif origin == "L":
                    j += 1
                    sequence1Str += gapChar
                    sequence2Str += str(self.sequence2[j])
                    midStr += indelChar
                elif origin == "D":
                    i += 1
                    j += 1
                    sequence1Str += str(self.sequence1[i])
                    sequence2Str += str(self.sequence2[j])
                    midStr += conservationChar if self.sequence1[i] == self.sequence2[j] else mutationChar
                
                indexHints1 += " "
                indexHints2 += " "
                if i > 0 and i % indexHintInterval == 0:
                    # Discard some characters at the end to make room for the index
                    # hint, and add the index hint
                    indexHints1 = indexHints1[:-len(str(i))] + str(i)
                    
                if j > 0 and j % indexHintInterval == 0:
                    indexHints2 = indexHints2[:-len(str(j))] + str(j)
            
            # Discard the first character, so that printed sequences start indexing
            # at 1 rather than 0 (because the index hints are shifted to the left)
            indexHints1 = indexHints1[1:]
            indexHints2 = indexHints2[1:]
            print("Solution No " + str(k + 1) + ":")
            
            for k in range(0, ((len(sequence1Str) - 1) // maxLineLength) + 1):
                indices = slice(k * maxLineLength, min((k + 1) * maxLineLength, len(sequence1Str)))
                print("           ", indexHints1[indices])
                print("Sequence 1:", sequence1Str[indices])
                print("           ", midStr[indices])
                print("Sequence 2:", sequence2Str[indices])
                print("           ", indexHints2[indices])
                print()
                
            print()

def main():
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

    aligner = Aligner(score, sequence1, sequence2, openPenalty, extendPenalty)

if __name__ == "__main__":
    cProfile.run("main()")
