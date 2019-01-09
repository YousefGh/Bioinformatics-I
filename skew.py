def skewCG(genome):
    '''
    skewness calculates the difference between 'C' and 'G' concentration
    This is connected to the place of the ori (where the cell replicates itself
    before cell division)

    :param genome: string sequence representing a DNA partition
    :return: an array indexing skewness of genome at each nucleotide
    '''
    skewness = [0] # i = 0
    skew = 0
    for i in genome:
        nucleotide = i.lower()
        if nucleotide == 'c':
            skew -= 1
        elif nucleotide == 'g':
            skew += 1
        skewness.append(skew)
    return skewness

def minSkewCG(genome):
    '''
    :param genome: string sequence representing a DNA partition
    :return: integers representing the minima of a skew array
    '''
    skewness = skewCG(genome)
    minimum = min(skewness)
    minima = []
    for index, nucleotide in enumerate(skewness):
        if nucleotide == minimum:
            minima.append(index)

    return minima