# Bioinformatics I 
#### Finding Hidden Messages in DNA (Bioinformatics I) by University of California San Diego

This course, [Bioinformatics I](https://www.coursera.org/learn/dna-analysis), focuses on the intersection between Computer Science and Biology (Bioinformatics). It takes the problem of finding the [Origin of Replication](https://en.wikipedia.org/wiki/Origin_of_replication) by applying different algorithms based on a biological concept that tells where the <i>ori</i> could be located, until we find a fairly fast algorithm.

### Files
Docs of each function is included in the source code, this is just an overview
* <i>cProfile_function.py</i> <br/>
  Detailed Calculation of the speed of a function execution

* <i>clumpFinder.py</i> <br/>
  Finds a pattern of length k (k-mer) in a genome sequence that is appearing many times (t) within a window L 

* <i>frequencyArray.py</i> <br/>
  An array that is mapping each occurring k-length pattern (k-mer) to array[patternToNumber] <br/>
  patternToNumber and numberToPattern is a way to convert nucleotide letters into numbers (Exactly as in number conversions into other bases where we have many bases like octal (base-8) decimal (base-10) binary (base-2)...). <br/>
  We have A-C-G-T represented as base-4 (Exactly similar to hexadecimal)

* <i>frequentPatterns.py</i> <br/>
  Counts the patterns of length k (k-mer) that are most occurring in a sequence

* <i>motifsFinder.py</i> <br/>
  Finds a pattern that is appearing in all sequences with the minimum number of mismatches
    so-called motifs. Details of each function are included in the source code
    
* <i>sequenceComplement.py</i> <br/>
  Finds the complement (reversed & complemented pattern) of a DNA string. DNA is usually a double-stranded, A is paired with T, and C is paired with G.
      
* <i>skew.py</i> <br/>
  skewness calculates the difference between 'C' and 'G' concentration in a DNA sequence. <i>ori</i> is located at where #C-#G is suddenly increasing after a decrease (It is not that easy to find this turning point)