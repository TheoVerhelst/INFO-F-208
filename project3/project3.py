from math import log
from copy import deepcopy
import os, sys
from enum import Enum

#sys.setrecursionlimit(1500)

### FILE project1.py ###

AminoAcid = Enum("AminoAcid", "A R N D C E Q G H I L K M F P S T W Y V B Z X")
# Make the str() function print only the amino acid letter by redefining __str__
AminoAcid.__str__ = lambda self: self.name[-1]

class Sequence:
    """Represents a sequence of amino acids.
    It is suited to compute the best alignment between two sequences."""
    
    def __init__(self, filename, sequence_id):
        """Constructor. Loads the sequence from a file.
        
        Parameters:
            filename: the filename of the file containing the sequence.
            sequence_id: the identifier of the sequence in the file.
                This is the string between the two first bars "|" on the line
                preceding the sequence.
            """
        self.load_from_file(filename, sequence_id)
    
    def __getitem__(self, index):
        """Returns the elements at the specified index. Allows the syntax
        sequence[i].
        """
        return self.sequence[index]
        
    def __len__(self):
        """Returns the lenght of the sequence. Allows the syntax len(sequence)."""
        return len(self.sequence)

    def load_from_file(self, filename, sequence_id):
        """Loads the sequence data from a file.
        
        Parameters:
            filename: the filename of the file containing the sequence.
            sequence_id: the identifier of the sequence in the file.
                This is the string between the two first bars "|" on the line
                preceding the sequence.
        """
        found_sequence = False
        sequence = ""
        with open(filename) as file:
            for line in file:
                if not found_sequence:
                    # We found the sequence if this is the right identifier line
                    found_sequence = (line[0] == ">"
                            and line.split("|")[1] == sequence_id)
                else:
                    if line[0] != ">":
                        sequence += line
                    else:
                        # Stop adding lines to the sequence
                        found_sequence = False
                        
        if sequence == "":
            raise RuntimeError("Sequence not found in the specified file.")

        try:
            sequence = [AminoAcid[character] for character in sequence                     if not character.isspace()]
        except KeyError as e:
            raise RuntimeError("An amino acid is not valid in the sequence: "
                    + str(e))

        # We didn't worked on self.sequence directly to preserve the internal
        # state in case of failure
        self.sequence = sequence
        
class Matrix:
    """Represents a m*n matrix.
    It is used to abstract the implementation of matrices, and allow possible
    optimisations easily, without having to change the whole code.
    All lines don't have to have the same length, so it may be used to represent
    a triangular matrix.
    
    """
    def __init__(self, m = 0, n = 0, value=None):
        """Construct the matrix
        
        Parameters:
            m, n: height and width of the matrix
            value: the initial value to put in the cells.
        """
        self.matrix = [[value] * n for i in range(m)]
    
    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self.matrix[key[0]][key[1]]
        else:
            return self.matrix[key]
            
    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            self.matrix[key[0]][key[1]] = value
        else:
            self.matrix[key] = value
    
    def __str__(self):
        cell_width = max([max([len(str(c)) for c in line]) for line in self.matrix])
        res = ""
        for line in self.matrix:
            for cell in line:
                res += str(cell).rjust(cell_width) + " "
            res += "\n"
        return res
            
    def append(self, iterable):
        """Adds a line to the bottom of the matrix."""
        self.matrix.append(list(iterable))
    
class Aligner:
    """Holds functions and variable related to the sequence alignment."""
    
    move_values = {"D" : (-1, -1), "T" : (-1, 0), "L" : (0, -1), "0" : (0, 0)}

    def __init__(self, sequence, profile,
            gap_penalties,
            sub_alignments = 0,
            display_all_solutions = True,
            max_line_length = 80):
        """Constructor, with all the parameters of the algorithm."""
        self.sequence = sequence
        self.profile = profile
        self.gap_penalties = gap_penalties
        self.sub_alignments = sub_alignments
        self.display_all_solutions = display_all_solutions
        self.solutions = []
        # The maximum number of characters displayed on a line
        self.max_line_length = max_line_length
        
        # Execute the algorithm
        self.fill_matrices()
        

    @staticmethod
    def make_move(coordinates, move_char, reverted = False):
        """Return the tuple coordinates, after that the specified movement has
        been applied to it.
        
        Parameters:
            coordinates: a tuple of integers
            move_char: one of the key of Aligner.move_values
            reverted: indicates whether the movement has to be negated before
                being applied.
        
        Return value: coordinates +- Aligned.move_values[move_char]
        """
            
        if not reverted:
            return coordinates[0] + Aligner.move_values[move_char][0],                     coordinates[1] + Aligner.move_values[move_char][1]
        else:
            return coordinates[0] - Aligner.move_values[move_char][0],                     coordinates[1] - Aligner.move_values[move_char][1]

    @staticmethod
    def find_maximum(matrix):
        """Finds s maximum in a matrix.
        
        Parameters:
            matrix: the matrix
        
        Return value: a tuple containing the (i, j) coordinates of the maximum
        """
        maximum = float("-inf")
        for i, line in enumerate(matrix):
            for j, value in enumerate(line):
                if value > maximum:
                    maximum = value
                    i_max, j_max = i, j
        return i_max, j_max
    
    def fill_matrices(self):
        """Fill the matrices in order to compute the alignment."""
        m, n = len(self.sequence) + 1, len(self.profile) + 1
        
        # Construct the alignment matrix first
        self.S = Matrix(m, n, 0)
        self.backtrace = Matrix(m, n, ["0"])
        
        # Compute the value of all cells of the matrix
        for i in range(1, m):
            for j in range(1, n):
                self.fill_cell(i, j)
        
        # The first cell of the solution is the maximum cell
        start_cell = self.find_maximum(self.S)
        self.find_alignment(*start_cell)
        
        #Print the result
        self.print_alignment()
        
        # Search for remaining subalignments
        for k in range(self.sub_alignments):
            self.clear_path(*start_cell, self.solutions[0]["path"])
            self.solutions = []
            start_cell = self.find_maximum(self.S)
            self.find_alignment(*start_cell)
            self.print_alignment(True, k)
        
    def fill_cell(self, i, j):
        """Computes the value in self.S and self.backtrace for the cell at
        position (i, j) according to the alignment rules.
        """
        
        # We use a dict to find the value and the origin string for the
        # current cell at the same time
        choices = {
            "T" : self.S[i - 1, j] + self.gap_penalties[j - 1],
            "L" : self.S[i, j - 1] + self.gap_penalties[j - 1],
            # We decrease the index for accessing the sequences because the
            # sequences have one less elements than the matrices
            "D" : self.S[i - 1, j - 1] + self.profile[j - 1, self.sequence[i - 1]],
            "0" : 0
        }
            
        # Find the maximum value in the dict
        self.S[i, j] = max(choices.values())
        
        # Add all origins that can lead to this maximum value in the backtrace
        self.backtrace[i, j] = ""
        for key in choices:
            if choices[key] == self.S[i, j]:
                self.backtrace[i,j] += key

    def clear_path(self, i, j, path):
        """Fills with zeros cells of the solution of a local alignment in order
        to compute the next subalignment. It also recomputes the cells impacted
        by the zeroing.
        
        Parameters:
            i, j: the coordinates of the first cell of the path
            path: the list of moves in this solution
        """
        cleared_cells = set()
        
        for move in path:
            self.S[i, j] = 0
            cleared_cells.add((i, j))
            i, j = self.make_move((i, j), move)
        first_cell = (i, j)
        
        for i in range(first_cell[0], len(self.sequence) + 1):
            for j in range(first_cell[1], len(self.profile) + 1):
                for origin in self.backtrace[i, j]:
                    if (i, j) not in cleared_cells                             and self.make_move((i, j), origin) in cleared_cells:
                        self.fill_cell(i, j)
                        cleared_cells.add((i, j))
        
        
    def find_alignment(self, i, j, current_solution = ""):
        """Recursively iterates on the backtrace matrix to find valid paths from
        the cell at position (i, j) to a cell with value 0.
        
        All paths meeting this criteria are appended to self.solutions (or just
        the first one if not self.display_all_solutions).
        
        Parameters:
            i, j: the coordinates of the currently explored cell
            current_solution: a string reprensenting the currently explored solution
        """
        # Stop if we are in a zero cell
        if "0" in self.backtrace[i, j]:
            self.solutions.append({"origin" : (i, j), "path" : current_solution})
        else:
            for possibility in self.backtrace[i, j]:
                self.find_alignment(*self.make_move((i, j), possibility),
                        current_solution + possibility)
                
                # Terminate after the first possiblity if we need only one solution
                if not self.display_all_solutions:
                    return
        
    def print_alignment(self, sub_alignment = False, sub_alignment_number = 0):
        """Prints all the alignments found so far in self.solutions.
        
        Parameters:
            sub_alignment: indicates whether we print a sub-alignment or not.
            sub_alignment_number: the index of the sub-alignment, if applicable.
        """
        gap_char = "-"
        indel_char = " "
        conservation_char = ":"
        mutation_char = "."
        # The interval between two index hints dipslayed above or below the sequence
        index_hint_interval = 10
        max_line_length = self.max_line_length - len("Sequence 1: ")
        
        if sub_alignment:
            print("Sub alignment No", str(sub_alignment_number + 1) + ":")
        else:
            print(len(self.solutions), "solution" + ("s" if len(self.solutions) > 1 else ""), "were found.")
            
        for k, solution in enumerate(self.solutions):
            sequence_str, profile_str, mid_str= "", "", ""
            index_hints_1, index_hints_2 = "", ""
            i, j = solution["origin"]
            i -= 1
            j -= 1
            
            # Iterate over the solution backward
            for origin in reversed(solution["path"]):
                i, j = self.make_move((i, j), origin, reverted = True)
                
                if origin == "T":
                    sequence_str += str(self.sequence[i])
                    profile_str += gap_char
                    mid_str += indel_char
                elif origin == "L":
                    sequence_str += gap_char
                    profile_str += str(self.profile.consensus[j])
                    mid_str += indel_char
                elif origin == "D":
                    sequence_str += str(self.sequence[i])
                    profile_str += str(self.profile.consensus[j])
                    mid_str += conservation_char if str(self.sequence[i])                             == self.profile.consensus[j] else mutation_char
                
                index_hints_1 += " "
                index_hints_2 += " "
                if i > 0 and i % index_hint_interval == 0 and origin != "L":
                    # Discard some characters at the end to make room for the index
                    # hint, and add the index hint
                    index_hints_1 = index_hints_1[:-len(str(i))] + str(i)
                    
                if j > 0 and j % index_hint_interval == 0 and origin != "T":
                    index_hints_2 = index_hints_2[:-len(str(j))] + str(j)
            
            # Discard the first character, so that printed sequences start indexing
            # at 1 rather than 0 (because the index hints are shifted to the left)
            index_hints_1 = index_hints_1[1:]
            index_hints_2 = index_hints_2[1:]
            
            if len(self.solutions) > 1:
                print("Solution No " + str(k + 1) + ":")
            
            for k in range(0, ((len(sequence_str) - 1) // max_line_length) + 1):
                indices = slice(k * max_line_length,
                        min((k + 1) * max_line_length,len(sequence_str)))
                print("          ", index_hints_1[indices])
                print("Sequence :", sequence_str[indices])
                print("          ", mid_str[indices])
                print("Consensus:", profile_str[indices])
                print("          ", index_hints_2[indices])
                print()
                
            print()

### END FILE project1.py ###

gap = "-"
aminoAcids = "ARNDCEQGHILKMFPSTWYV"

# The following probabilities has been gathered from the Swiss-prot database
# statistics page
p={ "E" : 0.0674, "H" : 0.0227, "S" : 0.066,  "K" : 0.0582, "Y" : 0.0292,
	"T" : 0.0535, "W" : 0.0109, "L" : 0.0965, "R" : 0.0553, "M" : 0.0241,
	"Q" : 0.0393, "A" : 0.0826, "V" : 0.0687, "F" : 0.0386, "C" : 0.0137,
	"N" : 0.0406, "P" : 0.0472, "I" : 0.0593, "D" : 0.0546, "G" : 0.0708}

logBase = 2

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
        
        self.computeConsensus()
    
    def __getitem__(self, t):
        return self.matrix[int(t[0])][str(t[1])]
    
    def __len__(self):
        return len(self.matrix)
    
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
    
    def computeConsensus(self):
        self.consensus = ""
        for column in self.matrix:
            self.consensus += max(column, key=lambda aa : column[aa])

def main():
    # Script parameters
    terminal_width = 80
    usedAligner = "muscle"
    sequence_file = "test.fasta"
    sequence_ids = ["D6C652", "P46935"]

    if usedAligner == "muscle":
        msaFilename = "msaresults-MUSCLE.fasta"
    else:
        msaFilename = "msaresults-CLUSTAL.fasta"
    
    pssm = Pssm(msaFilename)
    
    print(pssm)
    
    for sequence_id in sequence_ids:
        sequence = Sequence(sequence_file, sequence_id)
        print("*" * terminal_width)
        print("Alignment of", sequence_id)
        Aligner(sequence, pssm, gap_penalties = [-1] * len(pssm),
                display_all_solutions = False, sub_alignments = 3,
                max_line_length = terminal_width)
    
    print("PSSM consensus sequence:", pssm.consensus)

if __name__ == "__main__":
    main()
