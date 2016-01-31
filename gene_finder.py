# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: David Papp

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


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        The first four unit tests confirm that all of the letters work.
    	 The last one is to test that it returns 'Unknown' if it's some other letter.

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    >>>> get_complement('U')
    'Unknown'
    """

    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'
    else:
        return 'Unknown'



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        These given unit tests are sufficient because they contain all four types of nucleotides.
    
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    answer = ""
    for x in range(0, len(dna)):
        answer = answer + get_complement(dna[len(dna) - x - 1])
    return answer


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        The last unit test is to make sure that the function works when the start codon is immediately followed by a stop codon.
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGTAG")
    'ATG'

    """
    count = 0
    while count < len(dna):
    	current_frame = dna[count:count+3]
    	if current_frame == "TGA" or current_frame == "TAA" or current_frame == "TAG":
    		return dna[0:count]
    	count += 3

    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        The last unit test is to make sure that the function only starts reading from a start codon.

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TAAATGTAG")
    ['ATG']
    """
    answer = []
    loc = 0
    while loc < len(dna):
    	dna_substr = dna[loc:len(dna)]
    	if dna_substr[0:3] == "ATG":
    		new_sequence = rest_of_ORF(dna_substr)
    		answer.append(new_sequence)
    		loc += len(new_sequence)
    	else:
    		loc += 3
    return answer


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        This test is sufficient because it contains a start codon in each of the three possible frame locations.
        It also finds a stop codon in the third frame configuration, which confirms that the stop codons are also read correctly.
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    frame1 = find_all_ORFs_oneframe(dna)
    frame2 = find_all_ORFs_oneframe(dna[1:len(dna)])
    frame3 = find_all_ORFs_oneframe(dna[2:len(dna)])
    return frame1 + frame2 + frame3


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        This test is sufficient because it contains a start codon from both directions.
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    dna_reverse = get_reverse_complement(dna)
    normal = find_all_ORFs(dna)
    reverse = find_all_ORFs(dna_reverse)
    return normal + reverse


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this


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


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this

if __name__ == "__main__":
    import doctest
    #print find_all_ORFs_oneframe("ATGTAG")
    #print get_complement("A")
    #print get_reverse_complement("ATGCCCGCTTT")
    #print get_reverse_complement("CCGCGTTCA")
    #print rest_of_ORF("ATGTGAA")
    #print rest_of_ORF("ATGAGATAGG")
    #print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    #print find_all_ORFs("ATGCATGAATGTAG")
    #print find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals())
    #doctest.testmod()