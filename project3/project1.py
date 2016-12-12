from enum import Enum
from copy import deepcopy
import os

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
            sequence = [AminoAcid[character] for character in sequence \
                    if not character.isspace()]
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
        self.load_from_file(filename)
    
    def get_score(self, amino_acid_A, amino_acid_B):
        """Returns the match score between two amino acids according to
        the loaded matrix.
        """
        index_A = self.indices.index(amino_acid_A)
        index_B = self.indices.index(amino_acid_B)
        
        # The matrix is triangular, so we need to test the two possibilities for
        # indexing since one of the two would result in an index error
        if index_A > index_B:
            return self.matrix[index_A, index_B]
        else:
            return self.matrix[index_B, index_A]
    
    def load_from_file(self, filename):
        """Loads a substitution matrix from a file.
        
        Parameter:
            filename: the filename of the matrix to load. It must be a .iij file
            as given in the blosum archive.
        """
        self.matrix = Matrix()
        self.indices = []
        
        with open(filename) as file:
            found_amino_acid_list = False
            
            for line in file:
                # Skip commentary lines
                if line.strip().startswith("#"):
                    continue
                
                if not found_amino_acid_list:
                    # The first non-commentary line is the list of amino acids
                    # We reload this list because some files may order the acids
                    # differently
                    self.indices = [AminoAcid[character] for character in line.split()]
                    found_amino_acid_list = True
                else:
                    self.matrix.append([int(number) for number in line.split()])
    
class Aligner:
    """Holds functions and variable related to the sequence alignment."""
    
    move_values = {"D" : (-1, -1), "T" : (-1, 0), "L" : (0, -1), "0" : (0, 0)}

    def __init__(self, score, sequence_A, sequence_B,
            open_penalty = -12,
            extend_penalty = -2,
            local = True,
            sub_alignments = 0,
            display_all_solutions = True,
            max_line_length = 80):
        """Constructor, with all the parameters of the algorithm."""
        self.score = score
        self.sequence_A = sequence_A
        self.sequence_B = sequence_B
        self.open_penalty = open_penalty
        self.extend_penalty = extend_penalty
        self.local = local
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
            return coordinates[0] + Aligner.move_values[move_char][0], \
                    coordinates[1] + Aligner.move_values[move_char][1]
        else:
            return coordinates[0] - Aligner.move_values[move_char][0], \
                    coordinates[1] - Aligner.move_values[move_char][1]

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
        """Applies the Smith-Waterman or the Needleman-Wunsch algorithm,
        depending on wether we do local or global alignment.
        """
        m, n = len(self.sequence_A) + 1, len(self.sequence_B) + 1
        
        # Construct the alignment matrix first
        self.S = Matrix(m, n, 0)
        self.backtrace = Matrix(m, n, ["0"])
        
        if not self.local:
            for i in range(m):
                self.S[i, 0] = self.open_penalty + i * self.extend_penalty
                self.backtrace[i, 0] = ["T"]
                
            for j in range(1, n):
                self.S[0, j] = self.open_penalty + j * self.extend_penalty
                self.backtrace[0, j] = ["L"]
        
        self.V, self.W = deepcopy(self.S), deepcopy(self.S)
        
        # Compute the value of all cells of the matrix
        for i in range(1, m):
            for j in range(1, n):
                self.fill_cell(i, j)
        
        # The first cell of the solution is the maximum if we do local alignment,
        # or the bottom-right cell if we do global alignment
        start_cell = self.find_maximum(self.S) if self.local else (m - 1, n - 1)
        self.find_alignment(*start_cell)
        
        #Print the result
        self.print_alignment()
        
        # Search for remaining subalignments
        if self.local:
            for k in range(self.sub_alignments):
                self.clear_path(*start_cell, self.solutions[0]["path"])
                self.solutions = []
                start_cell = self.find_maximum(self.S)
                self.find_alignment(*start_cell)
                self.print_alignment(True, k)
        
    def fill_cell(self, i, j):
        """Computes the value in self.S, self.W, self.v and self.backtrace for
        the cell at position (i, j) according to the local or global alignment
        rules.
        """
        self.V[i, j] = max(self.S[i - 1, j] + self.open_penalty,
                self.V[i - 1, j] + self.extend_penalty)
        self.W[i, j] = max(self.S[i, j - 1] + self.open_penalty,
                self.W[i, j - 1] + self.extend_penalty)
        
        # We use a dict to find the value and the origin string for the
        # current cell at the same time
        choices = {
            "T" : self.V[i, j],
            "L" : self.W[i, j],
            # We decrease the index for accessing the sequences because the
            # sequences have one less elements than the matrices
            "D" : self.S[i - 1, j - 1] + self.score.get_score(
                    self.sequence_A[i - 1], self.sequence_B[j - 1])
        }
        
        if self.local:
            choices["0"] = 0
            self.V[i, j] = max(self.V[i, j], 0)
            self.W[i, j] = max(self.W[i, j], 0)
            
        # Find the maximum value in the dict
        self.S[i, j] = max(choices.values())
        
        # Add all origins that can lead to this maximum value in the backtrace
        self.backtrace[i, j] = "".join([key for key in choices \
                if choices[key] == self.S[i, j]])

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
            self.V[i, j] = 0
            self.W[i, j] = 0
            cleared_cells.add((i, j))
            i, j = self.make_move((i, j), move)
        first_cell = (i, j)
        
        for i in range(first_cell[0], len(self.sequence_A) + 1):
            for j in range(first_cell[1], len(self.sequence_B) + 1):
                for origin in self.backtrace[i, j]:
                    if (i, j) not in cleared_cells \
                            and self.make_move((i, j), origin) in cleared_cells:
                        self.fill_cell(i, j)
                        cleared_cells.add((i, j))
        
        
    def find_alignment(self, i, j, current_solution = ""):
        """Recursively iterates on the backtrace matrix to find valid paths from
        the cell at position (i, j) to a cell with value 0 if we do local
        alignement, or to the cell at position (0, 0) if we do global alignment.
        
        All paths meeting this criteria are appended to self.solutions (or just
        the first one if not self.display_all_solutions).
        
        Parameters:
            i, j: the coordinates of the currently explored cell
            current_solution: a string reprensenting the currently explored solution
        """
        
        # Stop if we are in a zero cell, or in the top-left cell
        if ("0" in self.backtrace[i, j] and self.local) \
                or (i == 0 and j == 0 and not self.local):
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
            sequence_A_str, sequence_B_str, mid_str= "", "", ""
            index_hints_1, index_hints_2 = "", ""
            i, j = solution["origin"]
            i -= 1
            j -= 1
            
            # Iterate over the solution backward
            for origin in reversed(solution["path"]):
                i, j = self.make_move((i, j), origin, reverted = True)
                
                if origin == "T":
                    sequence_A_str += str(self.sequence_A[i])
                    sequence_B_str += gap_char
                    mid_str += indel_char
                elif origin == "L":
                    sequence_A_str += gap_char
                    sequence_B_str += str(self.sequence_B[j])
                    mid_str += indel_char
                elif origin == "D":
                    sequence_A_str += str(self.sequence_A[i])
                    sequence_B_str += str(self.sequence_B[j])
                    mid_str += conservation_char if self.sequence_A[i] \
                            == self.sequence_B[j] else mutation_char
                
                index_hints_1 += " "
                index_hints_2 += " "
                if i > 0 and i % index_hint_interval == 0:
                    # Discard some characters at the end to make room for the index
                    # hint, and add the index hint
                    index_hints_1 = index_hints_1[:-len(str(i))] + str(i)
                    
                if j > 0 and j % index_hint_interval == 0:
                    index_hints_2 = index_hints_2[:-len(str(j))] + str(j)
            
            # Discard the first character, so that printed sequences start indexing
            # at 1 rather than 0 (because the index hints are shifted to the left)
            index_hints_1 = index_hints_1[1:]
            index_hints_2 = index_hints_2[1:]
            
            if len(self.solutions) > 1:
                print("Solution No " + str(k + 1) + ":")
            
            for k in range(0, ((len(sequence_A_str) - 1) // max_line_length) + 1):
                indices = slice(k * max_line_length,
                        min((k + 1) * max_line_length,len(sequence_A_str)))
                print("           ", index_hints_1[indices])
                print("Sequence 1:", sequence_A_str[indices])
                print("           ", mid_str[indices])
                print("Sequence 2:", sequence_B_str[indices])
                print("           ", index_hints_2[indices])
                print()
                
            print()

def main():
    # Script parameters
    scoring_matrix_filename = "blosum/blosum50.iij"
    sequence_A_file = "maguk-sequences.fasta"
    sequence_A_id = "O14936"
    sequence_B_file = "maguk-sequences.fasta"
    sequence_B_id = "Q86UL8"
    
    # Get the terminal width
    # If this line raises an error, just comment it and uncomment the following line
    terminal_width = int(os.popen('stty size', 'r').read().split()[1])
    #terminal_width = 80

    score = Score(scoring_matrix_filename)
    sequence_A = Sequence(sequence_A_file, sequence_A_id)
    sequence_B = Sequence(sequence_B_file, sequence_B_id)

    Aligner(score, sequence_A, sequence_B,
            open_penalty = -12,
            extend_penalty = -2,
            local = True,
            display_all_solutions = False,
            sub_alignments = 3,
            max_line_length = terminal_width)

if __name__ == "__main__":
    main()
