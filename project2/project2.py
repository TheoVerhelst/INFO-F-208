from os.path import splitext

def split_domain(domain_filename, blocks_number):
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

def main():
    # Generate PDZ-domain-A.fasta and PDZ-domain-B.fasta
    split_domain("PDZ-domain.fasta", 2)
    
    # Generate SH3-domain-A.fasta, SH3-domain-B.fasta, SH3-domain-C.fasta
    # and SH3-domain-D.fasta
    split_domain("SH3-domain.fasta", 4)

if __name__ == "__main__":
    main()
