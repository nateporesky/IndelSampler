# -*- coding: utf-8 -*-

# Test Change
from loguru import logger
import re
import sys
import os, os.path
from collections import Counter
from collections import defaultdict
from datetime import datetime
from math import floor
from statistics import median
import csv
import pandas as pd


# import pysam

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from io import StringIO



def main(params_dict, ampData, seqList, base_filename, IndelDF):


    #Extracting Amplicon Data
    seqRefUsed = params_dict['seqRefUsed']  # Index for the row we want to use
    Local_ID = ampData.iloc[seqRefUsed]["Identifier"]
    chr = ampData.iloc[seqRefUsed]["Chromosome"]
    Gene = ampData.iloc[seqRefUsed]["Gene"]
    Global_bp = ampData.iloc[seqRefUsed]["Position"]
    Ref_Seq = ampData.iloc[seqRefUsed]["Sequence"]



    #Extracting Parameters
    cutoff = params_dict['cutoff']
    DelLenCut = params_dict['DelLenCut']
    InsertMax = params_dict['InsertMax']
    Deletion = params_dict['Deletion']
    exclusion = params_dict['exclusion_threshold']
    gap_open = params_dict['nwalign_gap_open']
    gap_extend = params_dict['nwalign_gap_extend']
    
    param = {
        "nwalign_gap_open": gap_open,  # Currently -18, -40 was new, non-positive; default -20; -1 to -180 same alignment
        "nwalign_gap_extend": gap_extend,  # Currently -.6, -.6 was new, non-positive; default -0.6;  same alignment -19.9 to -0.001; less than gap_open
        "swalign_a": -8,
        "swalign_b": "nuc44m",
        "bamPattern": "shrimp2outARsorted.bam",
        "tile01startPos": 1,
        "tile02startPos": 1,
        "tile03startPos": 1,
    }



    # Scoring data file from ftp://ftp.ncbi.nih.gov/blast/matrices/
    # scoringdatafile=os.getenv('FLT3ITD_scoringMatrix')

    scoringdatafile = params_dict['scoringMatrix']

    seqs = seqList
    DelList = [] 
    IndexList = []
    current_date = datetime.today().strftime('%Y-%m-%d')
    logger.debug(f"Results: {current_date}\n")
    script_name = os.path.realpath(__file__)
    logger.debug(f"script: {script_name}")

    nuc44 = load_scoring_matrice(scoringdatafile)
    # print(nuc44)

    
    seqRef = Ref_Seq
    # Add in function to collect and then loop through more Reference sequences, currently hard coded in

    Cutoff = cutoff
    exclusion_threshold = exclusion  # Set exclusion threshold (e.g., exclude groups with less than 5% of total sequences)
    logger.info(f"Exclusion Threshold used: {exclusion_threshold * 100}%")

    for i in range(len(seqs)):  # Loops through all of the sample sequences, for each Reference sequence
        seqDel = seqs[i]
        firstIndex, DeletionSeq = AlignFunc(nuc44, seqRef, seqDel, DelList, Deletion, param, DelLenCut, InsertMax)
        IndexList.append(firstIndex)

    # significant_lengths, ConSeq, FinalConSeq, MSAscore, DeletionPoints = ConcensusDel(DelList, Cutoff, exclusion_threshold)  # Creates consensus seq
    significant_lengths, ConSeq, AverageIndexLength = ConcensusInDel(DelList, Cutoff, exclusion_threshold)  # Creates consensus seq


    VAF = round((len(DelList) / len(seqs)) * 100, 2)
    logger.info(f"Number of sequences containing duplications: {DupCounter}")

    logger.info("Consensus Sequences (Length, Count, Sequence):")
    for length, count in significant_lengths.items():
        logger.info(f"Length: {length}, Count: {count}, Insertion Sequence: {ConSeq[length][0]}, at Index: {round(ConSeq[length][1])}, QS = {ConSeq[length][2]}")


    logger.info(f"Reference Sequence Used: {Gene}, {chr}, {Global_bp}, {Ref_Seq}")
    logger.info(f"Sample Data Used: {filename}")
    logger.info(f"Number of Deletions/Insertions Found: {len(DelList)} from number of sequences: {len(seqs)}")
    logger.info(f"Average Index for each Length, starting from most common: {AverageIndexLength}")
    logger.info(f"VAF: {VAF}%")
    print(params_dict)
    # logger.info(f"Final Consensus Sequence with leading Insertion Point included: {FinalConSeq}")
    # logger.info(f"Avg Insertion Index from Clustered Sequences: {avg_first_index}")
    # logger.info(f"Final Insertion Point: {TrueIndex}")
    # logger.info(f"Final Score: {MSAscore}")

    #Use Pathlib
    output_folder = str(d_output)
    output_file = os.path.join(output_folder, f"{base_filename}_Ref-{Local_ID}_Out.txt")
    logger.info(output_file)
    if Deletion == True:
        status = "Deletion"
    else:
        status = "Insertion"



    with open(output_file, 'w') as f:
        f.write(f"File Ran: {os.path.basename(__file__)} as a {status} run\n")
        f.write(f"Results: {current_date}\n")
        f.write(f"Parameters used: Length Cutoffs: Min: {DelLenCut}, Max: {InsertMax}, Alignment Parameters: {param}\n")
        f.write(f"Exclusion Threshold used: {exclusion_threshold * 100}%\n")
        f.write(f"Reference Sequence Used: {Gene}, {chr}, {Global_bp}, {Ref_Seq}\n")
        f.write(f"Sample Data Used: {filename}\n")
        f.write(f"Number of Deletions/Insertions Found: {len(DelList)} from number of sequences: {len(seqs)}\n")
        f.write(f"VAF: {VAF}%\n")
        
        f.write(f"The Insertion/Deletion Sequences:\n")
        for length, count in significant_lengths.items():
            f.write(f"Length: {length-1}, Count: {count}, Insertion Sequence: {ConSeq[length][0][1:]}\nAt Index: {round(ConSeq[length][1])} at bp: {ConSeq[length][0][0]}, VAF = {round((count/len(seqs))*100, 3)}%, Quality Score: {round((ConSeq[length][2])*100, 3)}% similar.\n")
    
        # f.write(f"Final Consensus Sequence:\n")
        # f.write(f"Length: {len(FinalConSeq)-1}\nSequence: {FinalConSeq[1:]}\nInsertion Point: {avg_first_index} at {FinalConSeq[0]}\n")
        # f.write(f"Average Point of Insertion with respect to whole sequence: {TrueIndex}\n")
        # f.write(f"Final Score: {MSAscore}\n")

        output_file = os.path.join(output_folder, f"{base_filename}_Ref-{Local_ID}_Out.csv")
        logger.info(output_file)

        # create a new DataFrame with some data
        newIndelDF = pd.DataFrame({
        'Seqs': [filename],  # Make single values a list to match the expected structure
        'Date': [current_date],
        'Min_Length': [DelLenCut],
        'Max_Length': [InsertMax],
        "gap_open": [gap_open],
        "gap_extend": [gap_extend],
        "Exclusion_Threshild": [f"{exclusion_threshold * 100}%"],
        "Gene": [Gene],
        "Chromosome": [chr],
        "Global_Location": [Global_bp],
        "Reference_seq": [Ref_Seq],
        "seqFileName": [filename],
        "Num_InDels": [len(DelList)],
        "Num_seqs": [len(seqs)],
        "Con_Seq_Lengths": [list(significant_lengths.keys())],  # Convert dict keys to list
        "Con_Seq_Counts": [list(significant_lengths.values())],  # Convert dict values to list
        "VAF": [VAF],
        "Con_Seq": [str({key-1: value[0][1:] for key, value in ConSeq.items()})],
        "Con_Seq_Index": [str({key-1: round(value[1]) for key, value in ConSeq.items()})],  # Store as string or list if needed
        "Con_Seq_bp": [str({key-1: value[0][0] for key, value in ConSeq.items()})],
        "Con_Seq_VAF": [str({key-1: round((value / len(seqs)) * 100, 3) for key, value in significant_lengths.items()})],
        "Con_Seq_QS": [str({key-1: round(value[2] * 100, 3) for key, value in ConSeq.items()})],
    })

        return(newIndelDF)
    

def load_scoring_matrice(scoringdatafile):
    """
    Need to create scoring matrices for nwalign and for swalign
    
    """
    scorefile=open(scoringdatafile,'r')

    first=True
    print(f"scoringdatafile: {scoringdatafile}")
    for line in scorefile.readlines():
        #print(line)
        if line[0]=="#":
            continue
        if first:
            nuc44={}
            first=False
            columns=line.split()        
        else:
            rowdata=line.split()
            for vals in zip(rowdata[1:],columns):
                nuc44[(rowdata[0],vals[1])]=int(vals[0])
    
    return(nuc44)




'''
This uses swalign and then finds the largest inserted sequence
This finds the best cyclic permutation of the insertSeq against the refSeq, not currently used
From NP, to be worked on in future
- NP 8/22/24
'''

def findDupeSeq(insertSeq, refSeq, nuc44, gap_open, gap_extend):
    circPermLength=len(insertSeq)
    #maxSWscore=swalign(insertSeq, insertSeq)[0] #get a reference best score
    alignments =pairwise2.align.localds(insertSeq, insertSeq, nuc44, gap_open, gap_extend, one_alignment_only=True) #get a reference best score
    maxSWscore = alignments[0].score
    swScore=[]
    currInsertSeq=insertSeq
    # Go through various cyclic permutations of insert seq and compare to refseq
    for i in range(0,circPermLength):
        swScore.append(pairwise2.align.localds(currInsertSeq, refSeq, nuc44, gap_open, gap_extend, one_alignment_only=True).score)
        currInsertSeq=currInsertSeq[-1:]+currInsertSeq[:-1] # rotate insert sequence to the right

    # Find optimal permutation
    (bestScoreIndex,bestValue) = max(enumerate(swScore), key=(lambda x: x[1]))  
    dupeSeq = insertSeq[-bestScoreIndex:]+insertSeq[:-bestScoreIndex]
    
    # Find score and start vector
    [swScoreDupe,swStartDupe]=pairwise2.align.localds(dupeSeq, refSeq, nuc44, gap_open, gap_extend, one_alignment_only=True) 
    dupeLength=len(insertSeq)
    dupeStart=swStartDupe[1]-swStartDupe[0]+1+dupeLength
    dupeRelScore=float(swScoreDupe)/float(maxSWscore)
    return (dupeStart, dupeSeq, dupeRelScore, dupeLength)

# This uses swalign and then finds the largest deleted sequence

def findDeletionWithSW(seq1,seq2,nuc44, DelList, param, DelLenCut, InsertMax):

    #Request one of the best alignments using the scoring matrix nuc44, gapopen penalties of 20 and gap extend penalties of .5
    gap_open=param['nwalign_gap_open']
    gap_extend=param['nwalign_gap_extend']
     
    # https://biopython.org/docs/1.75/api/Bio.pairwise2.html
    one_align=True #Makes the output simpler
    alignments = pairwise2.align.localds(seq1, seq2, nuc44, gap_open, gap_extend, one_alignment_only=one_align)
    #print(f"parameters - gap_open: {gap_open}; gap_extend: {gap_extend}; one_aligment: {one_align} ")
    #alignments = pairwise2.align.globalms(seq1, seq2, matchScore, mismatchPenalty, gapOpen, gapExtend, one_alignment_only=one_align)
    #print(alignments)
    print()
    
    #Now we build a regex that will find all of the consecutive "-"'s
    p=re.compile('(-+)') 
    for i, a in enumerate(alignments):
        logger.debug('')
        logger.debug(f"alignment: {i}")
        logger.debug(f"seqA  = {a.seqA}")
        logger.debug(f"seqB  = {a.seqB}")
        logger.debug(f"score = {a.score}")
        logger.debug(f"start = {a.start}")
        logger.debug(f"end   = {a.end}")
        
        seqNuc = a.seqA # seqNuc refers to the reference sequence
        seqFind = a.seqB # seqFind refers to the sample sequence in which we expect the gap to occur
        del_find = re.findall(r'(-+)', seqFind) # Locates all strings of multiple gaps connected together

        # Added Catch to find seqs that don't have any deletions
        if not del_find:
            del_max = ''
            del_pos = 0
        else:
            del_max = max(del_find, key=len) # Reports the gap of longest length, the targeted deletion sequence
            del_pos = del_find.index(del_max) # Takes index of first gap location in gap span

        logger.debug(del_find)
        logger.debug(f"{del_pos}, {del_max}")
        
        # The p.finditer() call returns an iterator of all matches to the regex "p" in regSeqA.  
        # This is then used in a list comprehension to build up a list of tuples where the first element 
        # is the starting index and the 2nd is the length of the match.
        if (del_max == ''):
            del_list = [(0, 0)]
        else:
            del_list = [(m.start(),len(m.group())) for m in p.finditer(seqFind)] # Creates list of all deletions found in a sequence
        logger.debug(f"del_list: {del_list}")
        
        # This finds the largest element of the list of tuples, keyed on the 2nd tuple-element (i.e. length)
        (start, length) = max(del_list,key=(lambda x:x[1]))
        #print(f"start = {start}; length = {length}")
        indel_seq = seqNuc[start-1:start+length] # Forms the deletion sequence from starting point and length, and refers to reference sequence
        if len(indel_seq) > DelLenCut and len(indel_seq) < InsertMax: # Crops both too short and too long:
            DelList.append((start, indel_seq)) #Adds all of the sequenced deletions that are longer than the cut off to a list, cut off may vary depending on consistency of results
        
        #print(f"indel_seq = {indel_seq}")
    
    seq1A=alignments[0][0]
    seq2A=alignments[0][1]
    
    #print()
   # print(f"gap_open = {gap_open}")
    #print(f"gap_extend = {gap_extend}")

    print()
    #print(f"seq1A = {seq1A}")
    #print(f"seq2A = {seq2A}")

    # And then we iterate through the matching groups and find their length and location, and return the longest match
    # Includes catch for when no deletion is found
    if (del_max == ''):
        return (0, "")
    else:
        (firstIndex, length) = max([(m.start(), len(m.group())) for m in p.finditer(seq2A)], key=(lambda x: x[1])) # Searches for longest deletion in list
    #Now includes the base before and after the deletion sequence
  
    return (firstIndex, seq1A[firstIndex:firstIndex+length+1])

DupCounter = 0 # Global variable to count duplications
def findInsertsWithSW(readSeq,refSeq,nuc44, DelList, param, DelLenCut, InsertMax):
    '''
    Created from Deletion Script starting 8/8/2024
    - NP 8/22/24
    '''
    
    # Request one of the best alignments using the scoring matrix nuc44, gapopen penalties of 20 and gap extend penalties of .5
    gap_open=param['nwalign_gap_open']
    gap_extend=param['nwalign_gap_extend']
     
    # https://biopython.org/docs/1.75/api/Bio.pairwise2.html
    one_align=True #Makes the output simpler
    alignments = pairwise2.align.localds(readSeq, refSeq, nuc44, gap_open, gap_extend, one_alignment_only=one_align) #Same alignment program as deletions
    #print(f"parameters - gap_open: {gap_open}; gap_extend: {gap_extend}; one_aligment: {one_align} ")
    #alignments = pairwise2.align.globalms(seq1, seq2, matchScore, mismatchPenalty, gapOpen, gapExtend, one_alignment_only=one_align)
    print()
    
    #Now we build a regex that will find all of the consecutive "-"'s
    p=re.compile('(-+)') 
    for i, a in enumerate(alignments):
        logger.debug('')
        logger.debug(f"alignment: {i}")
        logger.debug(f"RefSeq  = {a.seqA}")
        logger.debug(f"ReadSeq  = {a.seqB}")
        logger.debug(f"score = {a.score}")
        
        seqNuc = a.seqA # seqNuc refers to the reference sequence and where we expect the gaps to occur
        seqFind = a.seqB # seqFind refers to the sample sequence, and where we expect the insert
        del_find = re.findall(r'(-+)', seqNuc) #Finds the gaps in Reference sequence in the indup version

        #Added Catch to find seqs that don't have any insertions
        if not del_find:
            del_max = ''
            del_pos = 0
        else:
            del_max = max(del_find, key=len)
            del_pos = del_find.index(del_max)

        logger.debug(del_find)
        logger.debug(f"{del_pos}, {del_max}")
        
        # The p.finditer() call returns an iterator of all matches to the regex "p" in regSeqA.  
        # This is then used in a list comprehension to build up a list of tuples where the first element 
        # is the starting index and the 2nd is the length of the match.
        if (del_max == ''): # Gives output in case of no insertions
            del_list = [(0, 0)]
        else:
            del_list = [(m.start(),len(m.group())) for m in p.finditer(seqNuc)] #Creates insertion list
        logger.debug(f"del_list: {del_list}")
        
        # This finds the largest element of the list of tuples, keyed on the 2nd tuple-element (i.e. length)
        (start, length) = max(del_list,key=(lambda x:x[1]))
        #print(f"start = {start}; length = {length}")
        
        indel_seq = seqFind[start-1:start+length]
        

        if len(indel_seq) > DelLenCut and len(indel_seq) < InsertMax: # Crops both too short and too long
            global DupCounter
            if seqFind.count(indel_seq) > 1: #Searches for exact duplications - NEED TO EDIT
                print("Duplication Found!")
                DupCounter += 1
            DelList.append((start, indel_seq)) #Adds all of the sequenced insertions that are longer than the cut off to a list

        #print(f"indel_seq = {indel_seq}")
    
    seq1A=alignments[0][0]
    seq2A=alignments[0][1]
    
    #print()
    #print(f"gap_open = {gap_open}")
    #print(f"gap_extend = {gap_extend}")

    print()
    #print(f"seq1A = {seq1A}")
    #print(f"seq2A = {seq2A}")


    # And then we iterate through the matching groups and find their length and location, and return the longest match
    ## Returns output in case of no Insertions
    if (del_max == ''):
        return (0, "")
    else:
        (firstIndex, length) = max([(m.start(), len(m.group())) for m in p.finditer(seq1A)], key=(lambda x: x[1])) 
    # Now includes the base before and after the deletion sequence
    return (firstIndex-1, seq2A[firstIndex:firstIndex+length+1])

def AlignFunc(nuc44, seqRef, seqDel, DelList, Deletion, param,  DelLenCut, InsertMax):
    '''
    Only slightly changed from original Align_kit function, for clarity purposes

    - NP 8/22/24
    '''

    funct_name = sys._getframe(0).f_code.co_name
    
    logger.debug("-------New Aligment--------")
    logger.debug(f"function: {funct_name}")
    
    print()
    logger.debug(f"seqRef = {seqRef}")
    logger.debug(f"seqInDel = {seqDel}")
    print()
    if Deletion == True:
        firstIndex, readSeqA = findDeletionWithSW(seqRef, seqDel, nuc44, DelList, param,  DelLenCut, InsertMax)
    else:
        firstIndex, readSeqA = findInsertsWithSW(seqRef, seqDel, nuc44, DelList, param,  DelLenCut, InsertMax)
    # Find the length of each sequence
    global Amp_Data
    global seqRefUsed
    indel_size = len(readSeqA)
        
    print()
    logger.debug(f"firstIndex = {firstIndex}")
    logger.debug(f"size = {indel_size}")
    logger.debug(f"readSeqA = {readSeqA}")
    print()
    return firstIndex, readSeqA

    
    

def fastaArray(filename, max_lines=15000): # can change max_lines to determine how many seqs you would like to process
    '''
    Creates a fasta formated array with all of the headers and their respective seuqences from a Fasta File

    - NP 7/20/24
    '''
    fasta = open(filename, 'r')
    entries = []  # empty list to store (header, sequence) tuples
    header = ''  # initialize header
    seq = ''  # initialize sequence
    line_count = 0  # track the number of lines read

    for line in fasta:
        line = line.strip()
        line_count += 1

        if line.startswith('>'):
            if header:  # check if header is not empty (not the first header)
                entries.append((header, seq))
            header = line[1:]
            seq = ''
        else:
            seq += line

        if line_count >= max_lines:
            break

    if header:
        entries.append((header, seq))

    return entries

def fastqArray(filename):
    '''
    Creates a fasta formated array with all of the headers and their respective sequences from any fastq-ish formatted file
    -No max_lines

    - NP 8/4/24
    '''
    valid_bases = {'A', 'T', 'C', 'G', 'N'}
    with open(filename, 'r') as fasta:
        entries = []  # List to store (header, sequence) tuples
        header = ''   # Initialize header
        seq = ''      # Initialize sequence
        expecting_sequence = False  # Flag to track if the next line should be a sequence

        for line in fasta:
            line = line.strip()  # Remove leading/trailing whitespace

            if line.startswith('@'):
                if header and seq:  # If there's an existing entry, save it
                    entries.append((header, seq))
                header = line[1:]  # Update header (remove '@')
                seq = ''           # Reset sequence for new entry
                expecting_sequence = True  # Expect a sequence line next

            elif expecting_sequence:
                # Process the sequence line if it's the immediate line after '@'
                if all(base in valid_bases for base in line):
                    seq = line  # Only take this line as the sequence
                else:
                    # If there are invalid characters, log a warning
                    print(f"Warning: Invalid characters found in sequence line: {line}")
                expecting_sequence = False  # Reset flag as we've processed the sequence

        # Append the last entry if any
        if header and seq:
            entries.append((header, seq))

    return entries

def getSeqList(fasta):
    '''
    Small mini function that pairs with fasta and fastq function lists
    - NP 7/20/24
    '''
    seqs = []
    for entry in fasta:
        seqs.append(entry[1])
        # Now 'sequence' contains the sequence corresponding to the current entry
    return seqs


def clusterLengths(numDataSort, lenDataSort): 
    '''
Function not changed from original, only made slight changes to integrate with rest of runction
- NP 7/24/24

clusterLengths function is designed to merge sequence lists that are within a difference of 1 from each other in length.

Initial Sorting:

The function first sorts numDataSort in descending order while keeping track of the original indices.
It then reorders lenDataSort to match the sorted order of numDataSort.
Merging Neighbors:

The function iterates twice (for immediate and secondary neighbors).
For each element in the sorted numDataSort, it checks its neighbors (both left and right) that are within a difference of 1 in lenDataSort.
If a neighbor is found, it merges the counts by adding the neighbor's count to the current element's count and sets the neighbor's count to 0.
Final Sorting:

After merging, the function sorts the lists again, this time in ascending order of lengths.
It reorders numDataSort to match the sorted order of lenDataSort.
Return Values:

The function returns the merged and sorted lists: numDataC (counts) and lenDataC (lengths).
Parameters
numDataSort (list of int): A list of counts corresponding to sequences.
lenDataSort (list of int): A list of lengths corresponding to sequences.
Returns
numDataC (list of int): The merged and sorted list of counts.
lenDataC (list of int): The sorted list of lengths.
This function is useful for clustering sequences that are close in length and aggregating their counts, which can help in reducing noise and improving the analysis of sequence data.
'''
    for i in [1, 2]:  # look for immediate and then secondary neighbors to merge
        (numSortIndex, numDataSort) = (list(t) for t in zip(*sorted(enumerate(numDataSort), key=lambda x: x[1], reverse=True)))
        lenDataSort = [lenDataSort[i] for i in numSortIndex]
        for iSort in range(0, len(numDataSort)):
            thisLen = lenDataSort[iSort]
            thisNum = numDataSort[iSort]
            if thisNum > 0:
                for j in [-i, i]:
                    # look to the left and right
                    leftLen = thisLen + j
                    leftIndex = next((x[0] for x in enumerate(lenDataSort) if x[1] == leftLen), None)
                    if leftIndex:
                        leftNum = numDataSort[leftIndex]
                        numDataSort[iSort] = numDataSort[iSort] + leftNum
                        numDataSort[leftIndex] = 0

    # sort in length order
    (lenDataCindex, lenDataC) = zip(*sorted(enumerate(lenDataSort), key=lambda x: x[1]))
    numDataC = [numDataSort[i] for i in lenDataCindex]

    return numDataC, lenDataC

'''
This is the main analysis function post alignment, ConcensusInDel
Creates a list of conesnsus sequences based on the lengths of the deletions/insertions in DelList
Created by NP, to be worked on in future
-NP 8/22/24

The ConcensusDel function processes a list of deletions to find consensus sequences based on their lengths.
It clusters the lengths, filters out outliers, aligns sequences of the same length, and generates a final consensus sequence.
If multiple consensus sequences are found, it performs multiple sequence alignment (MSA) to derive the final consensus.

Parameters
DelList (list of tuples): A list of deletions where each deletion is represented as a tuple (deletion point, sequence).
Cutoff (float): A cutoff percentage to select the top most frequent lengths.
exclusion_threshold (float): A threshold to exclude outlier groups based on their frequency.

Returns
significant_lengths (dict): A dictionary of significant lengths and their counts after exclusion.
consensus_sequences (dict): A dictionary of consensus sequences for each significant length.
final_consensus (str): The final consensus sequence obtained after MSA.
MSAscore (str): The score or status of the MSA process.
DeletionPoints (list): A list of deletion points corresponding to the sequences used.

Steps
Check for Empty Input:

If DelList is empty, the function returns an empty string.
Create Length Dictionary:

The function creates a dictionary LengthDict where the key is the length of the sequence and the value is a list of indices for each length.
Find Lengths and Count Occurrences:
It finds the length of each sequence in DelList and counts the occurrences of each length using Counter.

Cluster Lengths:
The function clusters the lengths using clusterLengths.
It updates length_counts with the clustered lengths.

Calculate Top Occurrences:
It calculates the number of top occurrences to display based on Cutoff.
It gets the top Cutoff% most frequent occurrences and logs the result.

Exclude Outlier Groups:
It excludes outlier groups based on exclusion_threshold and logs those lengths

Generate Consensus Sequences:
For each significant length, it filters sequences of the current valid length.
It aligns the sequences using pairwise2.align.globalms.
It builds a consensus sequence for each length by finding the most common base at each position.
Perform MSA if Needed:

If there are multiple consensus sequences, it performs MSA using multi_sequence_alignment.
If there is only one consensus sequence, it sets MSAscore to "No MSA Done".
'''
def ConcensusInDel(DelList, Cutoff, exclusion_threshold): 

    if not DelList: #Catches runs with no inserts
        return ""

    #MAKE DICTIONARY: {Key - Length : Value - [List of indices for each length]
    LengthDict = defaultdict(list)
    
    for i, item in enumerate(DelList):
        LengthDict[len(item[1])].append(i)

    # Find the length of each sequence
    lengths = [len(seq[1]) for seq in DelList]

    # Count the occurrences of each length
    length_counts = Counter(lengths)

    # Cluster the lengths
    numDataC, lenDataC = clusterLengths(list(length_counts.values()), list(length_counts.keys()))
    length_counts = dict(zip(lenDataC, numDataC))

    # Calculate the number of top occurrences to display
    Length_top = int(len(length_counts) * Cutoff)

    # Get the top Cutoff% most frequent occurrences
    Length_top_percent = dict(sorted(length_counts.items(), key=lambda item: item[1], reverse=True)[:Length_top])
    logger.debug(f"Top {Cutoff * 100}% lengths of all the Insertions: {Length_top_percent}")
    

    # Exclude outlier groups based on exclusion_threshold
    total_sequences = len(DelList)
    significant_lengths = {length: count for length, count in Length_top_percent.items() if (count / total_sequences) >= exclusion_threshold}
    logger.debug(f"Significant lengths after exclusion: {significant_lengths}")

    consensus_sequences = {}
    AverageIndexLength = []
    for length in significant_lengths:
        # Filter sequences of the current valid length
        # filtered_sequences = [seq[1] for seq in DelList if len(seq[1]) == length]
        filtered_sequences = list(map(lambda i: DelList[i][1], LengthDict[length]))
        DeletionPoints = list(map(lambda i: DelList[i][0], LengthDict[length]))

        if not filtered_sequences: 
            continue

        # Align sequences
        alignments = pairwise2.align.globalms(filtered_sequences[0], filtered_sequences[0], 2, -1, -0.5, -0.1, one_alignment_only=True)
        for seq in filtered_sequences[1:]:
            alignments = pairwise2.align.globalms(alignments[0][0], seq, 2, -1, -0.5, -0.1, one_alignment_only=True)

        # Build consensus sequence
        consensus_sequence = []
        quality_scores = []
        aligned_seq = alignments[0][0]

        for i in range(len(aligned_seq)):
            bases_at_position = [seq[i] for seq in filtered_sequences if i < len(seq)]
            base_counts = Counter(bases_at_position)
            if base_counts:
                consensus_base, count = base_counts.most_common(1)[0]
                # Calculate the proportion of bases that match the most common base
                match_proportion = count / len(bases_at_position)
                quality_scores.append(match_proportion)  # Store the quality score for this position

            else:
                continue
                #consensus_base = '-'  # Default to a gap if no bases are present
            consensus_sequence.append(consensus_base)

        quality_score = sum(quality_scores) / len(quality_scores)


        total = sum(DeletionPoints) / len(DeletionPoints)
        AverageIndexLength.append(total)


        consensus_sequences[length] = (''.join(consensus_sequence), total, quality_score)
        
        
    # if len(consensus_sequences) > 1:
    #     final_consensus, MSAscore, AverageIndexLength = multi_sequence_alignment(list(consensus_sequences.values()))
    # else:
        
    #     final_consensus = next(iter(consensus_sequences.values()))
    #     MSAscore = "No MSA Done"
    #     return significant_lengths, consensus_sequences, final_consensus, MSAscore, DeletionPoints
    
    return significant_lengths, consensus_sequences, AverageIndexLength

'''
Made 7/20/24
Not currently used as of 8/20/24
Was made for handling really rough sequences and alognments, not necessary currently after preprocessing
-NP
The multi_sequence_alignment function performs multiple sequence alignment (MSA) on a list of consensus sequences to find a final consensus sequence.
It uses ClustalW for the alignment process. 
Parameters
consensus_sequences (list of str): A list of consensus sequences to be aligned.
Returns
final_consensus (str): The final consensus sequence obtained after MSA.
stdout (str): The standard output from the ClustalW command.
Steps
Check for Empty Input:

If the input list consensus_sequences is empty, the function returns an empty string.
Write Sequences to a Temporary File:

The function writes the consensus sequences to a temporary file named consensus_sequences.fasta in FASTA format.
Run ClustalW:

The function constructs a ClustalW command line using ClustalwCommandline with the input file consensus_sequences.fasta.
It executes the ClustalW command and captures the standard output and error.
Parse the Aligned Sequences:

The function reads the aligned sequences from the file consensus_sequences.aln using AlignIO.read.
It initializes an empty list final_consensus to store the final consensus sequence.
Generate Final Consensus Sequence:

For each position in the alignment, the function counts the bases at that position across all sequences.
It determines the most common base at each position and appends it to final_consensus.
If no base is found, it appends a gap ('-').

Remove Leading and Trailing Gaps:
The function joins the list final_consensus into a string and removes leading and trailing gaps.

Clean Up Temporary Files:
The function removes the temporary files consensus_sequences.fasta, consensus_sequences.aln, and consensus_sequences.dnd.

Return Values:
The function returns the final consensus sequence and the standard output from the ClustalW command.
'''
def multi_sequence_alignment(consensus_sequences): # Uses MSA to find final concensus sequnece from all different lengths
    if not consensus_sequences:
        return ""

    # Write the consensus sequences to a temporary file
    with open("consensus_sequences.fasta", "w") as f:
        for i, seq in enumerate(consensus_sequences):
            f.write(f">seq{i}\n{seq}\n")
    

    # Run ClustalW

    # Execute the ClustalW command
    clustalw_cline = ClustalwCommandline("clustalw2", infile="consensus_sequences.fasta")
    stdout, stderr = clustalw_cline()  # Execute the ClustalW command

    # Parse the aligned sequences
    alignment = AlignIO.read("consensus_sequences.aln", "clustal")
    # Proceed with further analysis or processing of the alignment
    final_consensus = []

    for i in range(alignment.get_alignment_length()):
        bases_at_position = [record.seq[i] for record in alignment]
        base_counts = Counter(bases_at_position)
        if base_counts:
            consensus_base, _ = base_counts.most_common(1)[0]
        else:
            consensus_base = '-'
        final_consensus.append(consensus_base)

    # Remove leading and trailing gaps
    final_consensus = ''.join(final_consensus).strip('-')

    # Clean up temporary files
    os.remove("consensus_sequences.fasta")
    os.remove("consensus_sequences.aln")
    os.remove("consensus_sequences.dnd")

    return final_consensus, stdout


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="DEBUG")
    base_filename = "/base/"
    d_output = "/Outputs/"
    d_ampData = 'AmpData.xlsx'
    d_params = 'Params.xlsx'
    params_df = pd.read_excel(d_params)
    ampData = pd.read_excel(d_ampData)
    # Converting the parameters to a dictionary
    params_dict = pd.Series(params_df.Value.values, index=params_df.Parameter).to_dict()

    # Add a new handler that logs to stdout with level DEBUG
    # This is the known deletions seqs
    
    files = ["3'end.fasta", 'negative_control.fasta']
    
    #This functions extract sequences from FASTA/FASTQ files into a list, and for each file in the list executes main program
    IndelDF = pd.DataFrame()
    for file in files:
        filename = file
        logger.info(filename)
        if filename.lower().endswith('.fasta'):
            fasta = fastaArray(filename)
        elif filename.lower().endswith('.fastq'):
            fasta = fastqArray(filename)
        else:
            logger.error("This is not a valid file")
            
        seqs = getSeqList(fasta)

        SeqList = seqs

        newIndelDF = main(params_dict, ampData, SeqList, base_filename, IndelDF)
        IndelDF = pd.concat([IndelDF, newIndelDF])
        logger.info(newIndelDF)
        logger.info(IndelDF)
    IndelDF.to_csv(d_output + 'IndelDF.tsv', sep='\t',index=False)
    # pd.set_option('display.max_rows', None)  # Show all rows
    # pd.set_option('display.max_columns', None)  # Show all columns
    # print(IndelDF)  # or just `df` in a Jupyter notebook to display it


