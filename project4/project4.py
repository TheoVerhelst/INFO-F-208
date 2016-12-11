from enum import Enum

AminoAcid = Enum("AminoAcid", "A R N D C E Q G H I L K M F P S T W Y V B Z X")
# Make the str() function print only the amino acid letter by redefining __str__
AminoAcid.__str__ = lambda self: self.name[-1]

def parseCath(filename):
    res = []
    with open(filename) as file:
        for line in file:
            res.append({"file" : line[0:4] + ".dssp", "sequence" : line[4]})
            
    return res

def parseDssp(filename, sequence):
    res = []
    classes = \
    {
        "H" : "H",
        "G" : "H",
        "I" : "I",
        "E" : "E",
        "B" : "I",
        "T" : "T",
        "C" : "C",
        "S" : "C",
        ""  : "C"
    }
    
    with open(filename) as file:
        line = ""
        res = []
        
        while not line.strip().startswith("#"):
            line = file.readline()
            
        for line in file:
            identifier, aminoAcid, secondaryStruct = line.split()[2:5]
            if identifier.lower() == sequence.lower():
                # If there is no data for secondary structure, then the split
                # returned the next column, so we force to have empty string
                if secondaryStruct not in classes.keys():
                    secondaryStruct = ""
                res.append((AminoAcid[aminoAcid.upper()], classes[secondaryStruct.upper()]))
        return res

def writeSequence(filename, structure, identifier):
    with open(filename, "a") as file:
        file.write(">" + identifier + "|  |\n")
        aminoAcids = ""
        secondary = ""
        for aa, s in structure:
            aminoAcids += str(aa)
            secondary += s
        file.write(aminoAcids + "\n" + secondary + "\n")

dsspPrefix = "dataset/dssp_test/"
CathInfoFilename = "dataset/CATH_info_test.txt"
outputFilename = "out.fasta"

for parseData in parseCath(CathInfoFilename):
    structure = parseDssp(dsspPrefix + parseData["file"], parseData["sequence"])
