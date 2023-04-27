def read_file(file):
    file = read_file(file)
    FASTA_dict = {}
    FASTA_label = ""

    for line in file:
        if '>' in line:
            FASTA_label = line
            FASTA_dict[FASTA_label] = ""
        else:
            FASTA_dict[FASTA_label] += line

    print(FASTA_dict)
