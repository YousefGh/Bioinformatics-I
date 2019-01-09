def indexOfPattern(sequence, pattern):
    '''
    generates the starting indexes of each occurrence of pattern in a sequence

    :param sequence: Genome string
    :param pattern: k-mer
    :return: Array of indexes of pattern occurrence in a DNA sequence
    '''
    indexes = []
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i: i + len(pattern)] == pattern:
            indexes.append(i)

    return indexes



def patternToNumber(pattern):
    '''
    returns the number representing the pattern
    (clump finding problem helper function)

    :param num: a base 4 number
    :param k: length of the pattern to be converted
    :return: string pattern (k-mer)
    '''
    basesToNum = {'A': 0,'C': 1, 'G': 2,'T': 3} # "base" 4
    result = 0
    idx = len(pattern)  - 1

    for base in pattern:
        for dictBase, baseNum in basesToNum.items():
            if base == dictBase:
                result += pow(4, idx) * baseNum
        idx -= 1

    return result



def numberToPattern(num, k):
    '''
    the inverse of patternToNumber function. k = resulting string length
    (clump finding problem helper function)

    :param num: a base 4 number
    :param k: length of the pattern to be converted
    :return: string pattern (k-mer)
    '''
    basesToNum = {'A': 0,'C': 1, 'G': 2,'T': 3} # "base" 4
    result = ''
    for i in range(k):
        if num == 0:
            result = 'A' + result
        else:
            remainder = num % 4
            num = int(num / 4)
            for dictBase, baseNum in basesToNum.items():
                if baseNum == remainder:
                    result = dictBase + result
    return result




def frequencyArray(sequence, k):
    '''
    generates an array that is mapping each occurring k-length pattern to array(patternToNumber)

    :param sequence: DNA sequence string
    :param k: length of the pattern (k-mer)
    :return: Array of nucleotide bases as number indexes
    '''
    frequencyArray = [0] * 4**k

    for i in range(len(sequence) - k + 1):
        pattern = sequence[i:i + k]
        index = patternToNumber(pattern)
        frequencyArray[index] += 1

    return frequencyArray