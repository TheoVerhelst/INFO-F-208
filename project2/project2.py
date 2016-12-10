from os.path import splitext
from math import log2

amino_acids = "ARNDCQEGHILKMFPSTWYV"
# This dict speeds up the translation amino acid => index for amino_acids
# (at least in comparaison with amino_acid.index())
amino_acids_indices = dict((key, i) for i, key in enumerate(amino_acids))
gap_char = "X"
output_dir = "output/"

class CharMatrix:
    """Represents a triangular matrix indexed by a finite set of characters.
    It is suited to store values for relations between amino acids.
    """
    
    def __init__(self, keys, value):
        """Constructs the matrix.
        
        Parameters:
            keys: an iterable of characters, it will be the set of acceptable
                keys for indexing
            value: the initial value to put in the cells.
        """
        self.keys = "".join(list(keys))
        self.indices = dict((key, i) for i, key in enumerate(keys))
        self.matrix = [[value] * len(self.keys) for i in range(len(self.keys))]
    
    def __getitem__(self, key):
        return self.matrix[amino_acids_indices[key[0]]][amino_acids_indices[key[1]]]
            
    def __setitem__(self, key, value):
        self.matrix[self.indices[key[0]]][self.indices[key[1]]] = value
    
    def __str__(self):
        """Converts the matrix to a string."""
        
        cell_width = max(max(len(str(c)) for c in line) for line in self.matrix) + 1
        res = " ".join([c.rjust(cell_width) for c in self.keys]) + "\n"
        
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):
                if i >= j:
                    # We invert i and j in order to print a lower-half matrix
                    res += str(self.matrix[j][i]).rjust(cell_width) + " "
            res += "\n"
        return res
    
def split_domain(domain_filename, number_blocks):
    """Reads a .fasta file containing blocks of a domain, and writes the blocks
    in separate files. The blocks files have the same name as the domain file,
    except that a letter is appended just before the extension.
    
    The domain file must be formatted as explained in the statement of the
    project, thus copied from the BLOCKS database. Although the number of blocks
    could be detected from the first lines of the file, asking the user for this
    information make the code much clearer.
    
    Parameters:
        domain_filename: the filename of the domain file
        number_blocks: the number of differents blocks contained in the file
    
    Return value: the filenames of the created files
    """
    
    # Create the filenames of the block files
    root, extension = splitext(domain_filename)
    block_filenames = [output_dir + root + "-" + chr(ord("A") + i) + extension \
            for i in range(number_blocks)]
    
    # Open the block files
    block_files = [open(filename, "w") for filename in block_filenames]
    
    try:
        with open(domain_filename) as domain_file:
            for i, line in enumerate(domain_file):
                # local_index is 0 if the line is an identifier line,
                # or (1 + the current block) otherwise
                local_index = i % (number_blocks + 1)
                if local_index != 0:
                    block_files[local_index - 1].write(line)
    finally:
        for block_file in block_files:
            block_file.close()
    
    return block_filenames

def identity(string_a, string_b):
    """Computes the identity between two strings, in the range [0, 1]."""
    
    count = sum(1 if a == b else 0 for a, b in zip(string_a, string_b))
    return count / min(len(string_a), len(string_b))

def is_accepted_in_cluster(sequence, cluster, required_identity):
    """Checks if a sequence can be accepted in a cluster.
    Although the test is rather simple, we could think about more complicated
    criterion, this is why it is a separate function."""
    
    for other_sequence in cluster:
        if identity(other_sequence, sequence) >= required_identity:
            return True
    return False

def make_clusters(block_filename, required_identity):
    """Creates a list of clusters from the block file, and based on the given
    required identity. A cluster is represented by a list of strings.
    """
    
    clusters = []
    
    with open(block_filename) as file:
        for sequence in file:
            sequence = sequence.strip()
            found = False
            
            # Search for a cluster that can contains the current sequence
            for cluster in clusters:
                if is_accepted_in_cluster(sequence, cluster, required_identity):
                    # Add it if the cluster can accept it
                    cluster.append(sequence)
                    found = True
                    break
                    
            if not found:
                # Create a new cluster if there are no acceptable one
                clusters.append([sequence])
                
    return clusters
    
def weighted_frequencies(clusters):
    """Computes the weighted frequencies of all amino acids substitutions
    between the given clusters, i.e. the matrix f in the slides.
    """
    
    f = CharMatrix(amino_acids, 0)
    clusters_weights = [1 / len(cluster) for cluster in clusters]
    number_columns = len(clusters[0][0])
    number_clusters = len(clusters)
    counts = [[{} for j in range(number_columns)] for i in range(number_clusters)]
    
    for i, cluster in enumerate(clusters):
        for sequence in cluster:
            for c in range(number_columns):
                if sequence[c] != gap_char:
                    counts[i][c][sequence[c]] = counts[i][c].setdefault(
                            sequence[c], 0) + 1
    
    # Iterate over each possible pair of distinct clusters
    for i in range(number_clusters):
        for j in range(i + 1, number_clusters):
            weight = clusters_weights[i] * clusters_weights[j]
            # For each column of amino acid
            for col in range(number_columns):
                for a, frequency_a in counts[i][col].items():
                    for b, frequency_b in counts[j][col].items():
                        if amino_acids_indices[a] < amino_acids_indices[b]:
                            f[a, b] += weight * frequency_a * frequency_b
                        else:
                            f[b, a] += weight * frequency_a * frequency_b
    return f
    
def normalized_sum(matrices):
    """Computes the mean of the given matrices."""
    
    keys = matrices[0].keys
    res = CharMatrix(keys, 0)
    
    for a in keys:
        for b in keys[keys.index(a):]:
            res[a, b] = sum(matrix[a, b] for matrix in matrices)
    
    return res
    
def biological_probabilities(f):
    """Computes the probabilities of substitution in the biological pattern,
    i.e. the matrix q in the slides.
    """
    
    q = CharMatrix(amino_acids, 0)
    total_number_substitutions = 0
    
    for a in amino_acids:
        for b in amino_acids[amino_acids_indices[a]:]:
            total_number_substitutions += f[a, b]
    
    for a in amino_acids:
        for b in amino_acids[amino_acids_indices[a]:]:
            q[a, b] = f[a, b] / total_number_substitutions
    
    return q

def residue_probablities(q):
    """Computes the residue probabilities, i.e. the vector p in the slides."""
    
    p = [0] * len(amino_acids)
    
    for a in amino_acids:
        p[amino_acids_indices[a]] = q[a, a] + (sum(q[a, b] for b in amino_acids
                if a != b)) / 2
    
    return p
    
def random_probabilities(p):
    """Computes the alignment probabilities in the random pattern, i.e. the
    matrix e in the slides.
    """
    
    e = CharMatrix(amino_acids, 0)
    
    for a in amino_acids:
        for b in amino_acids[amino_acids_indices[a]:]:
            if a == b:
                e[a, b] = p[amino_acids_indices[a]] ** 2
            else:
                e[a, b] = p[amino_acids_indices[a]] * p[amino_acids_indices[b]] * 2

    return e

def log_odds_ratio(q, e):
    """Computes the log-odds ratio, i.e. the matrix s in the slides."""
    
    s = CharMatrix(amino_acids, 0)
    
    for a in amino_acids:
        for b in amino_acids[amino_acids_indices[a]:]:
            s[a, b] = 0 if q[a, b] == 0 else round(2 * log2(q[a, b] / e[a, b]))
    
    return s

def make_blosum(domain_filename, number_blocks, required_identity, print_info = False):
    """Constructs a BLOSUM matrix from a fasta file containing a domain.
    
    Parameters:
        domain_filename: the filename of the domain to use
        number_blocks: the number of BLOCKS in the domain
        required_identity: the identity required in a cluster
    
    Return value: the BLOSUM matrix
    """
    
    if print_info:
        print("Making BLOSUM" + str(round(required_identity * 100)), "from",
                domain_filename)
    filenames = split_domain(domain_filename, number_blocks)
    frequencies = []
    
    for filename in filenames:
        if print_info:
            print("Clustering", filename, "...")
        clusters = make_clusters(filename, required_identity)
        if print_info:
            print("Computing frequencies...")
        frequencies.append(weighted_frequencies(clusters))
        
    if print_info:
        print("Normalizing BLOCKS...")
    f = normalized_sum(frequencies)
    if print_info:
        print("Computing other matrices...")
    q = biological_probabilities(f)
    p = residue_probablities(q)
    e = random_probabilities(p)
    s = log_odds_ratio(q, e)
    
    return s
    
def main():
    print(make_blosum("PDZ-domain.fasta", 2, 0.4, True))
    print(make_blosum("SH3-domain.fasta", 4, 0.4, True))
    print(make_blosum("PDZ-domain.fasta", 2, 0.7, True))
    print(make_blosum("SH3-domain.fasta", 4, 0.7, True))

from time import clock
if __name__ == "__main__":
    t = clock()
    main()
    print(clock() - t)
