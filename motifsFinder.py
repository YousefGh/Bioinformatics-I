def hammingDistance(pattern_1, pattern_2):
    '''
    hamming distance is the difference between two patterns based on letters
    1 letter mismatch will add 1 to the hamming distance

    :param pattern_1: pattern in a genome
    :param pattern_2: pattern in a genome
    :return: integer representing hamming distance
    '''
    distance = 0
    for i in range(len(pattern_1)):
        if pattern_1[i] != pattern_2[i]:
            distance += 1
    return distance


def ApproximateHammingMatching(pattern, sequence, d):
    '''
    matches all similar patterns based on the hamming distance
    (letter differences)
    :param pattern: string k-mer
    :param sequence: long string (genome)
    :param d: number of differing letters
    :return: indexes array of all similar patters to pattern
    '''

    k = len(pattern)
    indexes = []
    for i in range(len(sequence) - k + 1):
        nextPattern = sequence[i: i + k]
        if hammingDistance(pattern, nextPattern) <= d:
            indexes.append(i)
    return indexes

def approximateFrequentPatterns(sequence, k, d):
    frequentPatterns = {}
    for i in range(len(sequence) - k + 1):
        pattern = sequence[i: i + k]
        for neighbor in neighbors(pattern,d):
            if neighbor not in frequentPatterns:
                frequentPatterns[neighbor] = 0
            frequentPatterns[neighbor] += 1

    maxFrequency = -float('inf')
    for pattern, frequency in frequentPatterns.items(): # TODO change to .values()
        if  frequency > maxFrequency:
            maxFrequency = frequency

    mostFrequentPatterns = []
    for pattern, frequency in frequentPatterns.items():
        if frequency == maxFrequency:
            mostFrequentPatterns.append(pattern)

    return mostFrequentPatterns


def neighbors(sequence, d):
    '''
    generates words that are differing than sequence by d letters called neighbors
    :param sequence: string of dna
    :param d: integer representing allowed distance of mismatches
    :return: array of neighbors (strings)
    '''

    if len(sequence) == 1: # Base case
        return ['A','C','G','T']

    result = []
    for neighbor in neighbors(sequence[1:],d):
        if hammingDistance(sequence[1:],neighbor) < d:
            result.extend(['A'+neighbor,'C'+neighbor,'G'+neighbor,'T'+neighbor])
        else:
            result.append(sequence[0] + neighbor)

    return result

def motifsEnumeration(dna, k, d):
    '''
    generates a list of strings that are not exactly appearing in all sequence of dna

    :param dna: an array of string each representing a sequence of DNA
    :param k: k-mer (length of the pattern)
    :param d: allowed distance of mismatches
    :return: array of motifs (k-mers that appear in all sequences of dna)
    '''

    motifs = set()
    for sequence in dna:
        for i in range(len(sequence) - k + 1):
            pattern = sequence[i: i + k]
            for neighbor in neighbors(pattern, d):
                if all(any(neighbor2 in seq for neighbor2 in neighbors(neighbor, d)) for seq in dna):
                    motifs.add(neighbor)

    return list(motifs)


def medianString(dna, k):
    '''
    generates the most frequent string of length k (k-mer) with
    minimum number of mismatches called median string

    :param dna: collection of DNA sequences
    :param k: length of the pattern (k-mer)
    :return: median string
    '''
    # setting up all possible k-mers
    firstKmer = ''
    for letter in range(k):
        firstKmer += 'A'
    kmers = neighbors(firstKmer, k)

    # finding median in motifs
    distance = float('inf')
    median = ''
    for kmer in kmers:
        localDistance = patternToDnaDistance(dna, kmer)
        if localDistance < distance:
            distance = localDistance
            median = kmer
    return median


def patternToDnaDistance(dna, pattern):
    '''
    calculates the total hamming distance between a pattern (k-mer) and
    each string in the DNA

    :param dna: collection of sequences of DNA
    :param pattern: k-mer
    :return: total hamming distance between pattern and all dna strings
    '''
    distance = 0
    k = len(pattern)
    for sequence in dna:
        localDistance = float('inf')
        for i in range(len(sequence) - k + 1):
            p = sequence[i: i + k]
            newDistance = hammingDistance(pattern, p)
            if newDistance < localDistance:
                localDistance = newDistance
        distance += localDistance
    return distance

def profileMostProbable(genome, k, profile):
    '''
    Takes a profile matrix to give the most probable k-mer to be
    generated based on the genome collection
    :param genome: list of DNA sequences
    :param k: length of the pattern k-mer
    :param profile: nucleotide probability matrix
    :return: Most probable k-mer
    '''
    nucleotides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    prob = -1
    mostProbPattern = 'None'
    for i in range(len(genome) - k + 1):
        pattern = genome[i: i + k]
        j = 0
        newProb = 1
        for nuc in pattern:
            newProb *= profile[j][nucleotides[nuc]]
            j += 1
        if newProb > prob:
            prob = newProb
            mostProbPattern = pattern
    return mostProbPattern


def profile(motifs, pseudoCounts = True):
    '''
    generates a matrix of probabilities where each row represents a
    nucleotide A-C-G-T respectively. Each column represents the
    index of the motif. row x column will be the probability of tha
    nucleotide (row) appearing at motif index (column)

    :param pseudoCounts: Laplace's Rule added (probability initialization
    starts from one instead of 0). True by default in this function unless changed

    :param motifs: list of DNA sequences
    :return: profile matrix
    '''
    motifs = [''.join(sequence) for sequence in zip(*motifs)]
    profileMatrix = []
    for column in motifs:
        nucleotides = ['A', 'C', 'G', 'T']
        row = []
        for nuc in nucleotides:
            # count the number of each nuc and divide it by the total
            if pseudoCounts:
                row.append(1.0 + (float(column.count(nuc)) / float(len(column))))
            else:
                row.append(float(column.count(nuc)) / float(len(column)))

        profileMatrix.append(row)

    return profileMatrix


def score(motifs):
    '''
    counts the number of the highest occurring nucleotide in each sequence of dna motifs
    and sums them up to give the score
    :param motifs: collection of DNA sequences
    :return: score of motifs
    '''
    motifs = [''.join(sequence) for sequence in zip(*motifs)]

    nucleotides = ['A', 'C', 'G', 'T']
    score = 0
    for column in motifs:
        nucCount = []
        for nuc in nucleotides:
            nucCount.append(column.count(nuc))
        colScore = len(column) - max(nucCount)
        score += (colScore)

    return score


def greedyMotifSearch(dna, k, t):
    bestScore = float('inf')
    bestMotifs = []

    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i: i + k]]

        for j in range(1, t):
            currentProfile = profile(motifs)
            mostProbKmer = profileMostProbable(dna[j], k, currentProfile)
            motifs.append(mostProbKmer)

        currentScore = score(motifs)
        if currentScore < bestScore:
            bestScore = currentScore
            bestMotifs = motifs

    return bestMotifs


if __name__ == '__main__':
    fileLinesArray = open('dataset_159_3.txt', 'r').read().split('\n')
    # #
    k = int(fileLinesArray[0])
    # dna = fileLinesArray[1].split(' ')


    dna = []
    for i in range(1, len(fileLinesArray)):
        dna.append(fileLinesArray[i])


    # genome = fileLinesArray[0]
    # k = int(fileLinesArray[1])
    #
    # profileMtrix = []
    # for i in range(2, len(fileLinesArray) ):
    #     row = fileLinesArray[i].split(' ')
    #     print(row)
    #     profileMtrix.append([float(row[j]) for j in range(k)])
    # # print(profileMtrix)
    #
    #
    # print(profileMostProbable(genome, k, profileMtrix))

    print(medianString(dna, k))