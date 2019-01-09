def seqComplement(pattern):
    '''
    finds the complement (reversed & complemented pattern) of a DNA string

    :param pattern: k-mer
    :return: complement of k-mer
    '''
    complement = ''
    complementaryBases = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for base in reversed(pattern):
        for dictBase, dictComp in complementaryBases.items():
            if base == dictBase:
                complement += dictComp

    return complement