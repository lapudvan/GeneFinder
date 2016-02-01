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
    nucleotide2 = nucleotide1.replace('A', 'P')
    nucleotide3 = nucleotide2.replace('T', 'A')
    nucleotide4 = nucleotide3.replace('P', 'T')
    nucleotide5 = nucleotide4.replace('C', 'S')
    nucleotide6 = nucleotide5.replace('G', 'C')
    nucleotide7 = nucleotide6.replace('S', 'G')
    return nucleotide7

    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    
    # this unit test is to test the other possiblities
    >>> get_complement('T')
    'A'
    # this unit test is to test the other possiblities
    >>> get_complement('G')
    'C'
    # this unit test is to test if it can handle many letters at once
    >>> get_complement('TAGC')
    'ATCG'
    """


def get_reverse_complement(dna):
    comp = get_complement(dna)
    reverse_comp = comp[::-1]
    return reverse_comp

    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    
    # I beleive this is substantial enough to check if it is taking the complement and reversing it because all letters are represented.
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """


def rest_of_ORF(dna):
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

    # This is to test the letter 'C' and the 'TAA' stop codon.
    >>> rest_of_ORF("ATGAGATCGTAAG")
    'ATGAGATCG'
    """

def find_ATG_Start(dna, startSearchingFrom):
        n = startSearchingFrom
        z = n + 3
        while z <= len(dna):
            if dna[n:z] == 'ATG':
                return n 
            else:
                n = n + 3
                z = z + 3
        return 'ERROR in find_ATG_START'

def find_all_ORFs_oneframe(dna):
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
    
    # This is to test a string that does not statr with a start codon.
    >>> find_all_ORFs_oneframe("TGGATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

def find_all_ORFs(dna):
    bList = []
    frame1 = find_all_ORFs_oneframe(dna[0:])
    bList = bList + frame1
    frame2 = find_all_ORFs_oneframe(dna[1:])
    bList = bList + frame2
    frame3 = find_all_ORFs_oneframe(dna[2:])
    bList = bList + frame3
    return bList

    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    
    # This tests if the function works with it starting on a not start codon.
    >>> find_all_ORFs("TTTATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """


def find_all_ORFs_both_strands(dna):
    strand1 = find_all_ORFs(dna)
    opposite = get_reverse_complement(dna)
    strand2 = find_all_ORFs(opposite)
    cList = []
    cList = cList + strand1
    cList = cList + strand2
    return cList

    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    
    #I wanted to test an example that has no result
    >>> find_all_ORFs_both_strands("AAAAAAAA")
    []
    """














 # This is the week 1 stopping point











def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()