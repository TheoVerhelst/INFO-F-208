from os.path import splitext

amino_acids = "ARNDCEQGHILKMFPSTWYV"

class CharMatrix:
    """Represents a triangular matrix indexed by a finite set of characters.
    It is suited to store values for relation between amino acids.
    """
    def __init__(self, keys, value):
        """Construct the matrix
        
        Parameters:
            keys: an iterable of characters, it will be the set of acceptable
                keys for indexing
            value: the initial value to put in the cells.
        """
        self.keys = "".join(list(keys))
        self.matrix = [[value] * len(self.keys) for i in range(len(self.keys))]
    
    def __getitem__(self, key):
        return self.matrix[self.keys.find(key[0])][self.keys.find(key[1])]
            
    def __setitem__(self, key, value):
        self.matrix[self.keys.find(key[0])][self.keys.find(key[1])] = value
    
    def __str__(self):
        """Converts the matrix to a string."""
        
        cell_width = max([max([len(str(c)) for c in line]) for line in self.matrix])
        res = (" " * cell_width) + " ".join([c for c in self.keys])
        
        for i, line in enumerate(self.matrix):
            res += str(self.keys[i]).rjust(cell_width)
            
            for cell in line:
                res += str(cell).rjust(cell_width) + " "
                
            res += "\n"
            
        return res
            
    def append(self, iterable):
        """Adds a line to the bottom of the matrix."""
        self.matrix.append(list(iterable))
    
def split_domain(domain_filename, blocks_number):
    """Reads a .fasta file containing blocks of a domain, and writes the blocks
    in separate files. The blocks files have the same name as the domain file,
    except that a letter is appended just before the extension.
    
    The domain file must be formatted as explained in the statement of the
    project, thus copied from the BLOCKS database. Although the number of blocks
    could be detected from the first lines of the file, asking the user for this
    information make the code much clearer.
    
    Parameters:
        domain_filename: the filename of the domain file
        blocks_number: the number of differents blocks contained in the file
    """
    
    # Create the filenames of the block files
    root, extension = splitext(domain_filename)
    block_filenames = [root + "-" + chr(ord("A") + i) + extension \
            for i in range(blocks_number)]
    
    # Open the block files
    block_files = [open(filename, "w") for filename in block_filenames]
    
    try:
        with open(domain_filename) as domain_file:
            for i, line in enumerate(domain_file):
                # local_index is 0 if the line is an identifier line,
                # or (1 + the current block) otherwise
                local_index = i % (blocks_number + 1)
                if local_index != 0:
                    block_files[local_index - 1].write(line)
    finally:
        for block_file in block_files:
            block_file.close()

def identity(string_a, string_b):
    """Computes the identity between two strings, in the range [0, 1]."""
    
    length = min(len(string_a), len(string_b))
    count = 0
    for i in range(length):
        if string_a[i] == string_b[i]:
            count += 1
    return count / length

def is_accepted_in_cluster(sequence, cluster, required_identity):
    """Check if a sequence can be accepted in a cluster."""
    
    identity = max([identify(sequence, sequence_2) for sequence_2 in cluster])
    return identity >= required_identity

def make_clusters(block_filename, required_identity):
    """Create a list of clusters from the block file, and based on the given
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
                    
            if not found:
                # Create a new cluster if there are no acceptable one
                clusters.append([sequence])
                
    return clusters
    
def weighted_frequencies(clusters):
    
    f = CharMatrix(amino_acids, 0)
    clusters_weights = [1 / len(cluster) for cluster in clusters]
    columns = len(clusters[0][0])
    
    # Iterate over each possible pair of amino acids
    for a in amino_acids:
        for b in amino_acids[amino_acids.find(a):]:
            f_ab = 0 # Work on a temp variable, to save matrix accesses
            
            # Iterate over each possible pair of distinct clusters
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    
                    # For each column of amino acid
                    for col in range(columns):
                        # Get all the amino acids of the current column in
                        # both clusters
                        col_i = "".join([seq[col] for seq in clusters[i]])
                        col_j = "".join([seq[col] for seq in clusters[j]])
                        
                        f_ab += clusters_weights[i] * clusters_weights[j] \
                                * col_i.count(a) * col_j.count(b)
                        f_ab += clusters_weights[i] * clusters_weights[j] \
                                * col_i.count(b) * col_j.count(a)
            f[a, b] = f_ab
    return f

def main():
    # Generate PDZ-domain-A.fasta and PDZ-domain-B.fasta
    split_domain("PDZ-domain.fasta", 2)
    
    # Generate SH3-domain-A.fasta, SH3-domain-B.fasta, SH3-domain-C.fasta
    # and SH3-domain-D.fasta
    split_domain("SH3-domain.fasta", 4)

if __name__ == "__main__":
    main()
