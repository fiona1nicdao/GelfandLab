import pandas as pd

def create_amino_acid_codon_table(protein_seq: str, dna_seq: str, output_file: str = None):
    """
    Creates a DataFrame matching amino acids to their corresponding DNA codons,
    including 3-letter amino acid codes.

    Parameters:
        protein_seq (str): Amino acid sequence (single-letter format)
        dna_seq (str): DNA sequence
        output_file (str, optional): If provided, saves the table as a CSV file

    Returns:
        pd.DataFrame
    """

    # 1-letter to 3-letter amino acid mapping
    aa_three_letter = {
        'A':'Ala', 'R':'Arg', 'N':'Asn', 'D':'Asp',
        'C':'Cys', 'E':'Glu', 'Q':'Gln', 'G':'Gly',
        'H':'His', 'I':'Ile', 'L':'Leu', 'K':'Lys',
        'M':'Met', 'F':'Phe', 'P':'Pro', 'S':'Ser',
        'T':'Thr', 'W':'Trp', 'Y':'Tyr', 'V':'Val',
        '*':'Stop'
    }

    # Clean inputs
    protein_seq = protein_seq.strip().upper()
    dna_seq = dna_seq.strip().upper()
    dna_seq = ''.join([base for base in dna_seq if base in "ATGC"])

    # Validate length
    expected_length = len(protein_seq) * 3
    if len(dna_seq) < expected_length:
        raise ValueError("DNA sequence is too short for the given protein sequence.")

    # Split DNA into codons
    codons = [dna_seq[i:i+3] for i in range(0, expected_length, 3)]

    # Convert to 3-letter amino acids
    three_letter_seq = [aa_three_letter.get(aa, 'Xaa') for aa in protein_seq]

    # Create DataFrame
    df = pd.DataFrame({
        "Position": range(1, len(protein_seq) + 1),
        "Amino Acid (1-letter)": list(protein_seq),
        "Amino Acid (3-letter)": three_letter_seq,
        "Codon": codons
    })

    # Save to CSV if requested
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"CSV file saved as {output_file}")

    return df
"""
def create_amino_acid_codon_table(protein_seq: str, dna_seq: str, output_file: str = None):

    Creates a DataFrame matching amino acids to their corresponding DNA codons.
    
    Parameters:
        protein_seq (str): Amino acid sequence (single-letter format)
        dna_seq (str): DNA sequence
        output_file (str, optional): If provided, saves the table as a CSV file
        
    Returns:
        pd.DataFrame


    # Clean inputs
    protein_seq = protein_seq.strip().upper()
    dna_seq = dna_seq.strip().upper()
    dna_seq = ''.join([base for base in dna_seq if base in "ATGC"])

    # Validate length
    expected_length = len(protein_seq) * 3
    if len(dna_seq) < expected_length:
        raise ValueError("DNA sequence is too short for the given protein sequence.")
    
    # Split DNA into codons
    codons = [dna_seq[i:i+3] for i in range(0, expected_length, 3)]

    # Create DataFrame
    df = pd.DataFrame({
        "Position": range(1, len(protein_seq) + 1),
        "Amino Acid": list(protein_seq),
        "Codon": codons
    })

    # Save to CSV if requested
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"CSV file saved as {output_file}")

    return df
"""
protien_pGEX_4T_1_UNC_76_MYC = "MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHNMLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALDVVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPKSDLVPRGSPEFMRDLGTKMAELKFEAPLAKFEETDEWGGCDFISSQNAINDTLNLNLKDSSAGGKPDATKLRLLEDAVRDAHVSKNGGGVVGGGAGSISPNCNTLQGGSLIEIGLSDVGLVPGEGAGVGLDGLEKRSLAVGDHVDNFTETFGGSLEDLVNTFDEKITKCFGNYEENVEELAPVQVRSQEEIMNECQMWWTITGNFGNILPIDWSKSYTRQMHMPTLNLGQNHTKQQQQNRNQQQQLHNQSHQAYPHTNGLGSGSGSGLDAQTPGDEFNDLTSEDEAVANDLDMHALILNGLNGDMDDQPIKTVEEVIKEIDDIMDEAESPLDEPETCDSEVIEKAREVLGAPLYAEKLQYLTTTQLNELYMEMEVLIQELSETLINELALRDELEFEKELKNSFISLLLAVQNKRRQYHVEKKRGKFQGPEPKYLTTVIPYHLENGTPNNQSLQVLIKILKAINEDSPTVPALLTDYILKVLCPTWRPLGGSGGGSGEQKLISEEDLAAAS"
dna_pGEX_4T_1_UNC_76_MYC = "atgtcccctatactaggttattggaaaattaagggccttgtgcaacccactcgacttcttttggaatatcttgaagaaaaatatgaagagcatttgtatgagcgcgatgaaggtgataaatggcgaaacaaaaagtttgaattgggtttggagtttcccaatcttccttattatattgatggtgatgttaaattaacacagtctatggccatcatacgttatatagctgacaagcacaacatgttgggtggttgtccaaaagagcgtgcagagatttcaatgcttgaaggagcggttttggatattagatacggtgtttcgagaattgcatatagtaaagactttgaaactctcaaagttgattttcttagcaagctacctgaaatgctgaaaatgttcgaagatcgtttatgtcataaaacatatttaaatggtgatcatgtaacccatcctgacttcatgttgtatgacgctcttgatgttgttttatacatggacccaatgtgcctggatgcgttcccaaaattagtttgttttaaaaaacgtattgaagctatcccacaaattgataagtacttgaaatccagcaagtatatagcatggcctttgcagggctggcaagccacgtttggtggtggcgaccatcctccaaaatcggatctggttccgcgtggatccccggaattcatgcgtgacctcggcaccaaaatggccgagctaaaattcgaggctccattggcaaagttcgaggaaacggacgagtggggcggctgcgacttcatctcgagtcagaatgccataaacgacacgctcaacctgaatctgaaggactcctccgccggcggcaagccggatgccaccaagctgcggctcctagaggatgctgtacgcgatgcgcacgtgagcaagaacggaggaggagtagtgggaggaggagctggctccattagccccaattgcaacacgctccagggtgggtcgcttattgaaatcggcctgtccgatgtgggattggtgccgggcgagggagccggcgttggtctggacggtctggagaagcgctccctggccgtcggcgaccatgtggacaacttcacggagaccttcggcggaagtctggaggatctggtcaatacgttcgacgaaaagataaccaagtgtttcggcaactacgaagagaacgtcgaggagctggctccggtgcaggtgcgcagccaggaggagatcatgaacgaatgccaaatgtggtggaccatcacaggaaactttggcaacattctgcccatcgactggtccaagtcgtacacccgtcagatgcacatgcccacactgaatctgggccagaaccacacgaagcagcagcagcagaatcgcaaccagcagcagcagttgcacaaccagagccatcaggcgtatccacacaccaacggtttgggatccggatccggatccggcctggatgcccagactccgggcgatgagttcaacgatctgaccagtgaggacgaggcggtggccaacgatttggacatgcacgccctcatccttaacggcttgaacggtgacatggacgaccagccgattaagacggtagaggaggtgatcaaggagatcgacgacatcatggacgaggccgagagcccgctggacgagccggagacctgcgattcggaggtgatcgagaaggcgcgagaggtgctgggcgcccccctctatgcggaaaaactgcaatatctgaccacgacacagctgaacgagctgtacatggagatggaggtgctgatccaggagctgagcgagacgttgatcaacgagctggccctgcgcgacgagctggagttcgagaaggagctgaagaactcgttcatctctctgctcctcgccgttcagaacaaaaggcggcagtaccatgtggagaagaagcgcggcaagttccagggtcccgagcccaagtacttgaccacagtcattccgtatcacctggagaacggcacgcccaacaatcagtctctgcaggtcctgatcaagattctcaaggccatcaacgaggacagcccgacggttccggccctccttacggactacatcctgaaggtgctctgccctacatggcggccgctgggtggttctggtggtggttctggtgaacaaaaactcatctcagaagaggatctggcggccgcatcg"

# create_amino_acid_codon_table(protien_pGEX_4T_1_UNC_76_MYC, dna_pGEX_4T_1_UNC_76_MYC, "unc76_myc.csv")

protein_pGEX_4T_1_KHCtail_HA = "MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHNMLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALDVVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPKSDLVPRGSPEFMREAENKKRTLEEQIDSLREECAKLKAAEHVSAVNAEEKQRAEELRSMFDSQMDELREAHTRQVSELRDEIAAKQHEMDEMKDVHQKLLLAHQQMTADYEKVRQEDAEKSSELQNIILTNERREQARKDLKGLEDTVAKELQTLHNLRKLFVQDLQQRIRKNVVNEESEEDGGSLAQKQKISFLENNLDQLTKVHKQLVRDNADLRCELPKLEKRLRCTMERVKALETALKEAKEGAMRDRKRYQYEVDRIKEAVRQKHLGRRGPQAQIAKPIRSGQGAIAIRGGGAVGGPSPLAQVNPVNSGRLGGSGGGSGYPYDVPDYAAAAS"
dna_pGEX_4T_1_KHCtail_HA = "atgtcccctatactaggttattggaaaattaagggccttgtgcaacccactcgacttcttttggaatatcttgaagaaaaatatgaagagcatttgtatgagcgcgatgaaggtgataaatggcgaaacaaaaagtttgaattgggtttggagtttcccaatcttccttattatattgatggtgatgttaaattaacacagtctatggccatcatacgttatatagctgacaagcacaacatgttgggtggttgtccaaaagagcgtgcagagatttcaatgcttgaaggagcggttttggatattagatacggtgtttcgagaattgcatatagtaaagactttgaaactctcaaagttgattttcttagcaagctacctgaaatgctgaaaatgttcgaagatcgtttatgtcataaaacatatttaaatggtgatcatgtaacccatcctgacttcatgttgtatgacgctcttgatgttgttttatacatggacccaatgtgcctggatgcgttcccaaaattagtttgttttaaaaaacgtattgaagctatcccacaaattgataagtacttgaaatccagcaagtatatagcatggcctttgcagggctggcaagccacgtttggtggtggcgaccatcctccaaaatcggatctggttccgcgtggatccccggaattcatgcgggaggcagagaacaagaagcgcacgctcgaggaacaaatcgattcgctgcgcgaggaatgcgccaagctcaaggccgccgagcacgtttccgccgttaacgccgaggagaaacagcgggctgaggagctgcgctccatgttcgattctcagatggacgagctacgcgaagcccacacccggcaggtgtccgagctccgggacgaaattgccgccaagcagcacgaaatggacgagatgaaggatgtccatcaaaagctgctcttggcgcaccaacagatgacggccgactacgagaaggtgcgccaggaggatgccgagaagtccagcgagcttcagaacatcatcctcaccaacgagcgtcgggagcaagcgcgcaaagacctcaagggcctggaggacacggtggccaaggagttgcagacgctacacaacctgcgaaaacttttcgttcaggatctacagcaacgaatccgaaagaatgtcgtaaacgaggagagcgaggaggacggtggatcactcgcgcagaaacagaagatttccttcttggagaacaacctcgaccagctgaccaaggtgcacaagcaattggtgcgggacaacgccgatctgcggtgcgagctgcccaagctggagaagcgtctacgctgtaccatggagcgggtgaaagctctggagacagcgctcaaggaggcgaaggagggcgcaatgcgggatcgcaagcgctaccaatacgaggtggaccgcatcaaggaggcggtgcgacagaagcatctgggcagacgtggcccacaggcacagatcgcaaagccgatccggtccggccaaggtgcaatcgccattcgtggtggtggcgccgttggaggaccatccccgctagcccaggttaatcctgtcaactcggGTCGACtgggtggttctggtggtggttctggtTACCCATACGATGTTCCAGATTACGCTgcggccgcatcg"

# create_amino_acid_codon_table(protein_pGEX_4T_1_KHCtail_HA, dna_pGEX_4T_1_KHCtail_HA, "kchtail_ha.csv")
protein_unc76_binding = "DNFTETFGGSLEDLVNTFDEKITKCFGNYEENVEELAPVQVRSQEEIMNECQMWWTITGNFGNILPIDWSKSYTRQMHMPTLNLGQNHTKQQQQNRNQQQQLHNQSHQAYPHTNGLGSGSGSGLDAQTPGDEFNDLTSEDEAVANDLDMHALILNGLNGDMDDQPIKTVEEVIKEIDDIMDEAESPLDEPETCDSEVIEKAREVLGAPLYAEKLQYLTTTQLNELYMEMEVLIQELSETLINELALRDELEFEKELKNSFISLLLAVQNKRRQYHVEKKRGKFQGPEPKYLTTVIPYHLENGTPNNQSLQVLIKILKAINEDSPTVPALLTDYILKVLCPTWR"
dna_unc76_binding = "gacaacttcacggagaccttcggcggaagtctggaggatctggtcaatacgttcgacgaaaagataaccaagtgtttcggcaactacgaagagaacgtcgaggagctggctccggtgcaggtgcgcagccaggaggagatcatgaacgaatgccaaatgtggtggaccatcacaggaaactttggcaacattctgcccatcgactggtccaagtcgtacacccgtcagatgcacatgcccacactgaatctgggccagaaccacacgaagcagcagcagcagaatcgcaaccagcagcagcagttgcacaaccagagccatcaggcgtatccacacaccaacggtttgggatccggatccggatccggcctggatgcccagactccgggcgatgagttcaacgatctgaccagtgaggacgaggcggtggccaacgatttggacatgcacgccctcatccttaacggcttgaacggtgacatggacgaccagccgattaagacggtagaggaggtgatcaaggagatcgacgacatcatggacgaggccgagagcccgctggacgagccggagacctgcgattcggaggtgatcgagaaggcgcgagaggtgctgggcgcccccctctatgcggaaaaactgcaatatctgaccacgacacagctgaacgagctgtacatggagatggaggtgctgatccaggagctgagcgagacgttgatcaacgagctggccctgcgcgacgagctggagttcgagaaggagctgaagaactcgttcatctctctgctcctcgccgttcagaacaaaaggcggcagtaccatgtggagaagaagcgcggcaagttccagggtcccgagcccaagtacttgaccacagtcattccgtatcacctggagaacggcacgcccaacaatcagtctctgcaggtcctgatcaagattctcaaggccatcaacgaggacagcccgacggttccggccctccttacggactacatcctgaaggtgctctgccctacatggcgg"

create_amino_acid_codon_table(protein_unc76_binding,dna_unc76_binding,"unc76bindingsite.csv")

