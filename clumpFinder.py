import cProfile

# calculates descriptively the time of a function execution
# add @do_cprofile above the function
def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func


# @do_cprofile
def clumpFinder(sequence, k, L, t):
    '''
    finds a pattern (k-mer) in a genome sequence that is appearing many times (t) within a window L
    slow don't use for large datasets (Use effClumpFinder instead)

    :param sequence: Genome
    :param k: length of the pattern (k-mer)
    :param L: Length of the window at which the pattern appears
    :param t: number of times the pattern appears at a window
    :return: array of clumps
    '''

    # Do not try to use my frequentPatterns() as it runs on O(n^2)
    print(k, L, t)
    frequentPatterns = []

    # 1. find all frequencies of all possible patterns

    # I used a hashed frequency array instead of the function frequencyArray()
    for i in range(len(sequence) - k + 1):
        pattern = sequence[i: i + k]

        if not any(fp['pattern'] == pattern for fp in frequentPatterns):
            frequentPatterns.append({'pattern': pattern, 'freq': 1, 'idx': [i]})
        else:
            for p in frequentPatterns:
                if p['pattern'] == pattern:
                    p['freq'] += 1
                    p['idx'] += [i]
    print('done hashing')

    # 2. remove any frequency < t
    for p in frequentPatterns[:]:
        if p['freq'] < t:
            frequentPatterns.remove(p)

    # 3. Add clumps (pattern repeated t times, within a window L)
    clumps = []
    for p in frequentPatterns:
        idx = 0
        found = False
        while idx < len(p['idx']) and not found:
            j = idx + 1
            numberOfComparedElements = []
            while j < len(p['idx']):
                window = p['idx'][j] - p['idx'][idx] - k
                currentFrequency = j - idx + 1
                if window <= L and currentFrequency >= t:
                    clumps.append(p['pattern'])
                    found = True
                    break
                j += 1
            idx += 1

    return clumps


# @do_cprofile
def chunkedClumpFinder(sequence, k, L, t):
    """"
    uses a set for duplicates, the dict is for frequencies only instead of an array of dicts
    separates each window length into an algorithm chunk (sliding window)
    3x faster than clumpFinder() but still slow (many redundant operations!)

    :param sequence: Genome
    :param k: length of the pattern (k-mer)
    :param L: Length of the window at which the pattern appears
    :param t: number of times the pattern appears at a window
    :return: set of clumps
    """

    frequentPatterns = set([])
    for i in range(len(sequence)):
        window = sequence[i:i + L]
        frequencies = {}

        for j in range(len(window)):
            pattern = window[j:j + k]
            if pattern not in frequencies:
                frequencies[pattern] = 1
            else:
                frequencies[pattern] += 1
        for p in frequencies:
            if frequencies[p] >= t:
                frequentPatterns.add(p)
    return frequentPatterns


def fastClumpFinder(sequence, k, L, t):
    """"
    uses a set for duplicates, the dict is for frequencies only instead of an array of dicts
    separates each window length into an algorithm chunk (sliding window)
    3x faster than clumpFinder() but still slow (many redundant operations!)

    :param sequence: Genome
    :param k: length of the pattern (k-mer)
    :param L: Length of the window at which the pattern appears
    :param t: number of times the pattern appears at a window
    :return: set of clumps
    """

    # to be implemented ;)
    pass