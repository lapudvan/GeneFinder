# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Lauren Pudvan

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide1):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T') # this unit test is to test the other possiblities
    'A'
    >>> get_complement('G') # this unit test is to test the other possiblities
    'C'
    >>> get_complement('TAGC') # this unit test is to test if it can handle many letters at once
    'ATCG'
    """
    nucleotide2 = nucleotide1.replace('A', 'P')
    nucleotide3 = nucleotide2.replace('T', 'A')
    nucleotide4 = nucleotide3.replace('P', 'T')
    nucleotide5 = nucleotide4.replace('C', 'S')
    nucleotide6 = nucleotide5.replace('G', 'C')
    nucleotide7 = nucleotide6.replace('S', 'G')
    return nucleotide7


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT") # I beleive this is substantial enough to check if it is taking the complement and reversing it because all letters are represented.
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    comp = get_complement(dna)
    reverse_comp = comp[::-1]
    return reverse_comp


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATCGTAAG") # This is to test the letter 'C' and the 'TAA' stop codon.
    'ATGAGATCG'
    """
    length = len(dna)
    n = 0
    z = 3
    while z < length:
        if dna[n:z] == 'TAG' or dna[n:z] == 'TGA' or dna[n:z] == 'TAA':
            return dna[:n]
        else:
            n = n + 3
            z = z + 3
    return dna


def find_ATG_Start(dna, startSearchingFrom):
        """ Finds the first ATG in frame. This was made to assist the function find_all_ORFs_oneframe.
        """
        n = startSearchingFrom
        z = n + 3                       #frames of 3
        while z <= len(dna):
            if dna[n:z] == 'ATG':
                return n 
            else:
                n = n + 3
                z = z + 3
        return 'ERROR in find_ATG_START'


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TGGATGCATGAATGTAGATAGATGTGCCC") # This is to test a string that does not statr with a start codon.
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    aList = []
    location = 0
    while location <= len(dna)+3:
        h = find_ATG_Start(dna, location)
        if h != 'ERROR in find_ATG_START':
            j = rest_of_ORF(dna[h:])
            aList.append(j)
            end = h + len(j)
            location = end + 3
        else:
            return aList


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATGTAG']
    >>> find_all_ORFs("TTTATGCATGAATGTAG") # This tests if the function works with it starting on a not start codon.
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATGTAG']
    """
    bList = []
    frame1 = find_all_ORFs_oneframe(dna[0:])
    bList = bList + frame1
    frame2 = find_all_ORFs_oneframe(dna[1:])
    bList = bList + frame2
    frame3 = find_all_ORFs_oneframe(dna[2:])
    bList = bList + frame3
    return bList


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("AAAAAAAA") #I wanted to test an example that has no result
    []
    """
    strand1 = find_all_ORFs(dna)
    opposite = get_reverse_complement(dna)
    strand2 = find_all_ORFs(opposite)
    cList = []
    cList = cList + strand1
    cList = cList + strand2
    return cList


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("CCCATGCGAATGTAGCATCAAA") # More possiblities
    'ATGCTACATTCGCATGGG'
    >>> longest_ORF("CCCATGCGAATGTAGCATCAAAGTA") # More possiblities
    'ATGCTACATTCGCATGGG'
    """
    allORFList = find_all_ORFs_both_strands(dna)
    return max(allORFList, key=len)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 

        This can not be tested with doctests because of the random results. 
        To test this I ran it with one trial and had it print the randomDNA and the ListOfLongORFs.
        Then I did the calcualtions by hand to see if it was returning the length of the longest ORF.
        """
    ListOfLongORFs = []
    for i in range(num_trials):
        randomDNA = shuffle_string(dna)
        longString = longest_ORF(randomDNA)
        ListOfLongORFs.append(longString)
    return len(max(ListOfLongORFs, key=len))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGCCCGCTTTA") # Testing for a string with a multiple of 3 letters
        'MPAL'
        >>> coding_strand_to_AA("CCCGCTTTA") # Testing to see if doesn't start with ATG
        'PAL'
    """

    length = len(dna)
    n = 0
    stringOfAA = ''
    while (n + 3) <= length:
        tripplet = dna[n:n+3]
        amino_acid = aa_table[tripplet]
        stringOfAA  += amino_acid
        n = n + 3
    return stringOfAA


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        Cannot use doctests because the threshhold is not predictable. 
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFs = find_all_ORFs_both_strands(dna)
    listOfORFsLongEnough = []
    listOfAA = []
    for ORF in ORFs:
        if len(ORF) >= threshold:
            listOfORFsLongEnough.append(ORF)
    for longORF in listOfORFsLongEnough:
        AA = coding_strand_to_AA(longORF)
        listOfAA.append(AA)
    return listOfAA




if __name__ == "__main__":
    import doctest
    doctest.testmod()
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print gene_finder(dna)