import cProfile

# calculates descriptively the time of function execution
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

def frequencyCounter(sequence, pattern):
    '''
    counts the number of times the pattern occurs in a sequence

    :param sequence: genome
    :param pattern: k-mer
    :return: the frequency of pattern sequence in DNA
    '''
    counter = 0
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i: i + len(pattern)] == pattern:
            counter+=1

    return counter


def frequentPatterns(sequence, k):
    '''
    counts the patterns of length k (k-mer) that are most occurring in a sequence

    :param sequence: genome
    :param k: length of pattern (k-mer)
    :return: the most k length frequent patterns, k-mer "motif"
    '''
    frequentPatterns = []
    highestFreq = 0

    for i in range(len(sequence) - k):
        pattern = sequence[i:i + k]
        frequency = frequencyCounter(sequence, pattern)

        # if freq is higher and it should not be already added
        if frequency >= highestFreq and not any(fp['pattern'] == pattern for fp in frequentPatterns):
            highestFreq = frequency
            frequentPatterns.append({'pattern': pattern, 'freq': frequency})
            for p in frequentPatterns[:]:
                if p['freq'] < highestFreq:
                    frequentPatterns.remove(p)

    return frequentPatterns