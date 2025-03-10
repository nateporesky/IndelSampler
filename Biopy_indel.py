# -*- coding: utf-8 -*-
"""
source /opt/venvs/envOFA2/bin/activate

envOFA2) (base) wcrowe@urmc-sh.rochester.edu@smdmdlvm01:/mnt/ClinNGS_Share/Software_dev/InDelDuper/0.2.10/indelduper$ python3 scripts/biopython/biopy_KIT_del.py

Test

"""
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
from pathlib import Path


# import pysam

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from io import StringIO


def main(files, d_output, d_scoreMatrix, d_ampData): 
    '''
    Parameters:
    files: A list of files to process. Each file can be a string (for TSV files) or a tuple (for FASTA/FASTQ files).
    d_output: The output directory to save the output files.
    d_scoreMatrix: The path to the scoring matrix file.
    d_ampData: The path to the amplicon data file.
    '''

    #Extracting Amplicon Data
    ampExcel = pd.read_excel(d_ampData)
    IndelDF = pd.DataFrame()
    for file in files:
        logger.info(file)
        if isinstance(file, str) and file.lower().endswith('.tsv'): #This is the currently used workflow
            #AmpID_List = tsvArray(file)
            Amp_DF = pd.read_csv(file, sep='\t')
            # Group by the categorical variable
            Amp_DF = Amp_DF.groupby('AmpID')
            # Loop through each group
            for AmpID, seqs in Amp_DF:

                matching_rows = ampExcel[ampExcel['amp_id'] == AmpID].index
                AmpID_DF = Amp_DF.get_group(AmpID)
                for rownum in matching_rows:
                    ampData = ampExcel.iloc[[rownum]]
                    amp_dict = {}    

                    amp_dict['scoreMatrix'] = d_scoreMatrix #MIGHT NEED TO CHANGE
                    amp_dict['chr'] = ampData["chr"].iloc[0] #Used in output
                    amp_dict["gene"] = ampData["gene"].iloc[0] #Used in output
                    amp_dict["amp_id"] = ampData["amp_id"].iloc[0] #Used in output
                    amp_dict["hg19s"] = ampData["hg19s"].iloc[0] #Used in output
                    amp_dict["ref_seq"] = ampData["ref_seq"].iloc[0] #Used in output and then later alignment



                    #Extracting Parameters
                    amp_dict['cutoff'] = ampData['cutoff'].iloc[0] #Used in alignment
                    amp_dict['DelLenCut'] = ampData['DelLenCut'].iloc[0] #Used in alignment
                    amp_dict['InsertMax'] = ampData['InsertMax'].iloc[0] #Used in alignment
                    amp_dict['Deletion'] = bool(ampData['Deletion'].iloc[0]) #Used in alignment
                    amp_dict['exclusion'] = ampData['exclusion'].iloc[0] #Used in ConsensusInDel, Set exclusion threshold (e.g., exclude groups with less than 5% of total sequences)
                    amp_dict['gap_open'] = ampData['gap_open'].iloc[0] #Used in alignment
                    amp_dict['gap_extend'] = ampData['gap_extend'].iloc[0] #Used in alignment
                    # Retrieve the deletion status for the current row
                    deletion_status = ampExcel.loc[rownum, 'Deletion']
                    
                    # Set the text based on the deletion status
                    status_text = "del" if deletion_status else "ins"
                    
                    # Incorporate the status into the filename
                    filename = f"{AmpID}_{status_text}_{Path(file).name}" #Filename is made from ampID, indel type of run, and the original filename
                    
                    output_sum = d_output / f"{filename}_indelRunSum.txt" #filename used to save unique summary output
                    output_full = d_output / f"{filename}_indelRunFull.txt" #filename used to save unique full output

                    with open(output_sum, 'w') as f:
                        f.write(f"File Ran: {filename}\n")
                    
                    with open(output_full, 'w') as f:
                        f.write(f"File Ran: {filename}\n")

                    newIndelDF = findIndel(amp_dict, AmpID_DF, filename, output_sum, output_full)

                    IndelDF = pd.concat([IndelDF, newIndelDF])
                    logger.info(newIndelDF)
                    logger.info(IndelDF)
            IndelDF.to_csv(d_output / f'{Path(file).name}_IndelDF.tsv', sep='\t',index=False) #Save the complete IndelDF to a .tsv file


        elif isinstance(file, tuple) and file[0].lower().endswith('.fasta'): #This is the old workflow
            filename = file[0]
            rownum = file[1]
            fasta = fastaArray(filename)
            
            ampData = ampExcel.iloc[[rownum]]   
            seqs = getSeqList(fasta)
            filename = Path(filename).name
            SeqList = seqs
            output_sum = d_output / f"{filename}_indelRunSum.txt"
            output_full = d_output / f"{filename}_indelRunFull.txt"
            
            with open(output_sum, 'w') as f:
                f.write(f"File Ran: {filename}\n")
            
            with open(output_full, 'w') as f:
                f.write(f"File Ran: {filename}\n")

            newIndelDF = findIndel(ampData, SeqList, filename, output_sum, output_full)

            IndelDF = pd.concat([IndelDF, newIndelDF])
            logger.info(newIndelDF)
            logger.info(IndelDF)
            IndelDF.to_csv(d_output / 'IndelDF.tsv', sep='\t',index=False)

        elif isinstance(file, tuple) and file[0].lower().endswith('.fastq'):
            filename = file[0]
            amp_id = file[1]
            fasta = fastqArray(filename)
            for row in ampExcel.iterrows():
                if row[1]['amp_id'] == amp_id:
                    rownum = row[0]
                    ampData = ampExcel.iloc[[rownum]] # type: ignore
                        
                    seqs = getSeqList(fasta)
                    filename = Path(filename).name
                    SeqList = seqs
                    output_sum = d_output / f"{filename}_indelRunSum.txt"
                    output_full = d_output / f"{filename}_indelRunFull.txt"

                    with open(output_sum, 'w') as f:
                        f.write(f"File Ran: {filename}\n")
                    
                    with open(output_full, 'w') as f:
                        f.write(f"File Ran: {filename}\n")

                    newIndelDF = findIndel(ampData, SeqList, filename, output_sum, output_full)

                    IndelDF = pd.concat([IndelDF, newIndelDF])
                    logger.info(newIndelDF)
                    logger.info(IndelDF)
            IndelDF.to_csv(d_output / 'IndelDF.tsv', sep='\t',index=False)
        else:
            logger.error("This is not a valid file")
        
    #This functions extract sequences from TSV/FASTA/FASTQ files into a list, and for each file in the list executes main program
    # pd.set_option('display.max_rows', None)  # Show all rows
    # pd.set_option('display.max_columns', None)  # Show all columns
    # print(IndelDF)  # or just `df` in a Jupyter notebook to display it


def findIndel(amp_dict, AmpID_DF, filename, output_sum, output_full):
    '''
    Takes Parameter DF data, list of seqs, directory name, name of sample sequence file, and both initialized outputfiles

    '''
    Deletion = amp_dict['Deletion']
    chr = amp_dict['chr'] #Used in output
    gene = amp_dict["gene"] #Used in output
    amp_id = amp_dict["amp_id"] #Used in output
    hg19s = amp_dict["hg19s"] #Used in output

    status = "Deletion" if Deletion is True else "Insertion" #Used in to tell Alignfunc which to run and output
    print(status)
    print(type(Deletion))

    #Change for address of NUC.4.4
    scoringdatafile = amp_dict['scoreMatrix']
    
    DelList = [] #VERY IMPORTANT, this is where all InDels are stored
    IndexList = [] #VERY IMPORTANT, this is where all Indexes are stored for each InDel, indexes do not account for gap spaces being present
    current_date = datetime.today().strftime('%Y-%m-%d')
    logger.info(f"Results: {current_date}\n")

    nuc44 = load_scoring_matrice(scoringdatafile)
    # print(nuc44)

    
    seqRef = amp_dict["ref_seq"] # Format used in alignment
    Cutoff = amp_dict['cutoff']
    DelLenCut = amp_dict['DelLenCut'] #Used in alignment
    InsertMax = amp_dict['InsertMax'] #Used in alignment
    exclusion = amp_dict['exclusion'] #Used in ConsensusInDel, Set exclusion threshold (e.g., exclude groups with less than 5% of total sequences)
    gap_open = amp_dict['gap_open'] #Used in alignment
    gap_extend = amp_dict['gap_extend'] #Used in alignment

    # Add in function to collect and then loop through more Reference sequences, currently hard coded in
    logger.info(f"Exclusion Threshold used: {exclusion * 100}%")
    for i in range(len(AmpID_DF)):  # Loops through all of the sample sequences, for each Reference sequence
        seqInDel = AmpID_DF.loc[i, 'seq']
        seqHeader = AmpID_DF.loc[i, 'head']
        firstIndex, DeletionSeq = AlignFunc(nuc44, seqRef, seqInDel, seqHeader, DelList, Deletion, gap_open, gap_extend, DelLenCut, InsertMax, output_full)
        IndexList.append(firstIndex)

    # significant_lengths, ConSeq, FinalConSeq, MSAscore, DeletionPoints = ConcensusDel(DelList, Cutoff, exclusion_threshold)  # Creates consensus seq
    significant_lengths, ConSeq, AverageIndexLength = ConcensusInDel(DelList, Cutoff, exclusion)  # Creates consensus seq
    if significant_lengths == 0:
        newIndelDF = pd.DataFrame({
                'Seqs': filename,  # Make single values a list to match the expected structure
                'Date': current_date,
                'Min_Length': DelLenCut,
                'Max_Length': InsertMax,
                "gap_open": gap_open,
                "gap_extend": gap_extend,
                "Exclusion_Threshild": f"{exclusion * 100}%",
                "Gene": gene,
                "Amp_ID": amp_id,
                "Chromosome": chr,
                "Global_Location": hg19s,
                "Reference_seq": seqRef,
                "Num_InDels": [0],
                "Num_seqs": len(AmpID_DF),
                "Con_Seq_Lengths": [0],  # The current length
                "Con_Seq_Counts": [0],  # The count of this specific length
                "VAF": [0], # VAF for all InDels Detected
                "Con_Seq": [0],  # The consensus sequence without the first character
                "Con_Seq_Index": [0],  # The insertion index for the consensus sequence
                "Con_Seq_bp": [0],  # The bp (first character) for this consensus sequence
                "Con_Seq_VAF": [0],  # VAF for the consensus
                "Con_Seq_QS": [0],  # Quality Score for the consensus sequence
            })
        with open(output_sum, 'a') as f:
            f.write("No Significant Sequences found\n")
        return newIndelDF

    VAF = round((len(DelList) / len(AmpID_DF)) * 100, 2)
    

    logger.info("Consensus Sequences (Length, Count, Sequence):")
    for length, count in significant_lengths.items():
        logger.info(f"Length: {length}, Count: {count}, Insertion Sequence: {ConSeq[length][0]}, at Index: {round(ConSeq[length][1])}, QS = {ConSeq[length][2]}")


    # Calculates the average starting index value from only the sequences used to calculate the consensus sequences
    #TrueIndex = round(int(Amp_Data[seqRefUsed][2]) + avg_first_index)

    # logger.info(f"Reference Sequence Used: {Gene}, {chr}, {Global_bp}, {Ref_Seq}")
    # logger.info(f"Sample Data Used: {filename}")
    # logger.info(f"Number of Deletions/Insertions Found: {len(DelList)} from number of sequences: {len(seqs)}")
    # logger.info(f"Average Index for each Length, starting from most common: {AverageIndexLength}")
    # logger.info(f"VAF: {VAF}%")
    # logger.info(f"Final Consensus Sequence with leading Insertion Point included: {FinalConSeq}")
    # logger.info(f"Avg Insertion Index from Clustered Sequences: {avg_first_index}")
    # logger.info(f"Final Insertion Point: {TrueIndex}")
    # logger.info(f"Final Score: {MSAscore}")

    #Use Pathlib
    
    logger.info(output_sum)


    with open(output_sum, 'a') as f: #Output Summary File
        f.write(f"{status} run\n")
        f.write(f"Results: {current_date}\n")
        f.write(f"Parameters used: Length Cutoffs: Min: {DelLenCut}, Max: {InsertMax}, Gap Open: {gap_open}, Gap Extend: {gap_extend}\n")
        f.write(f"Exclusion Threshold used: {exclusion * 100}%\n")
        f.write(f"Reference Sequence Used: {gene}, {chr}, {hg19s}, {seqRef}\n")
        f.write(f"Sample Data Used: {filename}\n")
        f.write(f"Number of Deletions/Insertions Found: {len(DelList)} from number of sequences: {len(AmpID_DF)}\n")
        f.write(f"VAF: {VAF}%\n")
        
        f.write(f"The InDel Sequences:\n")
        for length, count in significant_lengths.items():
            f.write(f"Length: {length-1}, Count: {count}, InDel Sequence: {ConSeq[length][0][1:]}\nAt Index: {round(ConSeq[length][1])} at bp: {ConSeq[length][0][0]}, VAF = {round((count/len(AmpID_DF))*100, 3)}%, Quality Score: {round((ConSeq[length][2])*100, 3)}% similar.\n")
            f.write("\n")
        for length, count in significant_lengths.items(): #This grabs example alignments for each of the significant lengths
            target_seq = str(ConSeq[length][0][1:]).strip()  # Define target string
            target = f"InDel = {target_seq}"
            f.write(f"Target: {target}\n")
            target_found = False  # Flag to track if target was found

            # Open output_full for reading
            with open(output_full, 'r') as f2:
                buffer = []  # Initialize a buffer to store lines

                for line in f2:
                    # Add the current line to the buffer
                    buffer.append(line)

                    # Check if the current line contains the target
                    if target in line:
                        target_found = True  # Set flag to True if target is found
                        # Write the last 13 lines (12 previous + target line) to output_sum
                        f.writelines(buffer[-16:])  # Get the last 16 lines
                        break  

            # If target was not found, write a message to output_sum
            if not target_found:
                f.write(f"No exact match was located for target: {target}\n")
        # f.write(f"Final Consensus Sequence:\n")
        # f.write(f"Length: {len(FinalConSeq)-1}\nSequence: {FinalConSeq[1:]}\nInsertion Point: {avg_first_index} at {FinalConSeq[0]}\n")
        # f.write(f"Average Point of Insertion with respect to whole sequence: {TrueIndex}\n")
        # f.write(f"Final Score: {MSAscore}\n")


        # create a new DataFrame with some data
        newIndelDF = pd.DataFrame()
        # Loop through each ConSeq item (each represents a unique consensus sequence)
        for key, value in ConSeq.items():
            # Create a row for each Con_Seq
            row = {
                'Sample': filename,  # Make single values a list to match the expected structure
                'Date': current_date,
                'Del/Ins': status,
                'Min_Length': DelLenCut,
                'Max_Length': InsertMax,
                "gap_open": gap_open,
                "gap_extend": gap_extend,
                "Exclusion_Threshild": f"{exclusion * 100}%",
                "Gene": gene,
                "Amp_ID": amp_id,
                "Chromosome": chr,
                "Global_Location": hg19s,
                "Reference_seq": seqRef,
                "Num_InDels": len(DelList),
                "Num_seqs": len(AmpID_DF),
                "Con_Seq_Lengths": key-1,  # The current length
                "Con_Seq_Counts": significant_lengths[key],  # The count of this specific length
                "VAF": VAF,
                "Con_Seq": value[0][1:],  # The consensus sequence without the first character
                "Con_Seq_Index": round(value[1]),  # The insertion index for the consensus sequence
                "Con_Seq_bp": value[0][0],  # The bp (first character) for this consensus sequence
                "Con_Seq_VAF": round((significant_lengths[key]) / len(AmpID_DF) * 100, 3),  # VAF for the consensus
                "Con_Seq_QS": round(value[2] * 100, 3),  # Quality Score for the consensus sequence
            }
            
            # Append the row to the DataFrame
            newIndelDF = pd.concat([newIndelDF, pd.DataFrame([row])], ignore_index=True)

        return newIndelDF


def AlignFunc(nuc44, seqRef, seqInDel, seqHeader, DelList, Deletion, gap_open, gap_extend,  DelLenCut, InsertMax, output_full):
    '''
    Only slightly changed from original Align_kit function, for clarity purposes
    This is the parent function that calls either the deletion or the insertion function, and writes to full log output file
    - NP 8/22/24
    '''
    
    # https://www.oreilly.com/library/view/python-cookbook/0596001673/ch14s08.html
    funct_name = sys._getframe(0).f_code.co_name

    with open(output_full, 'a') as f:
        f.write(f"Deletion?: {Deletion} run\n")
        f.write("-------New Alignment--------\n")
        f.write(f"Sequence Header: {seqHeader}\n")
        f.write("\n")
        f.write(f"Parameters used: Length Cutoffs: Min: {DelLenCut}, Max: {InsertMax}, Gap Open: {gap_open}, Gap Extend: {gap_extend}\n")
        f.write("\n")
        f.write(f"seqRef   = {seqRef}\n")
        f.write(f"seqInDel = {seqInDel}\n")
        f.write("\n")

    if Deletion is True:
        firstIndex, readSeqA = findDeletionWithSW(seqRef, seqInDel, nuc44, DelList, gap_open, gap_extend, DelLenCut, InsertMax, output_full)
    else:
        firstIndex, readSeqA = findInsertsWithSW(seqRef, seqInDel, nuc44, DelList, gap_open, gap_extend, DelLenCut, InsertMax, output_full)
    # Find the length of each sequence
    indel_size = len(readSeqA[:-1])
        
    with open(output_full, 'a') as f:
        f.write(f"firstIndex = {firstIndex}\n")
        f.write(f"size = {indel_size}\n")
        f.write(f"InDel = {readSeqA[:-1]}\n")
        f.write("\n")
        
    return firstIndex, readSeqA

def load_scoring_matrice(scoringdatafile): #Untouched by NP, loads in Nuc.4.4 scoring matrix
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


def findDeletionWithSW(seq1,seq2,nuc44, DelList, gap_open, gap_extend, DelLenCut, InsertMax, output_full):
    '''
    Modified and integrated from WEC Script starting 7/6/2024

    Input: Reference and Sample sequence, The DelList to store found Indels, Parameters, and full outputfile

    Output: The first index of the deletion is returned, the sequence of the deletion in DelList, logged in output_full
    - NP 8/22/24
    '''
   
    one_align=True #Makes the output simpler
    alignments = pairwise2.align.globalds(seq1, seq2, nuc44, gap_open, gap_extend, one_alignment_only=one_align) #This creates the alignment, DOES NOT STORE INDEL, that is done with DelList
    #print(f"parameters - gap_open: {gap_open}; gap_extend: {gap_extend}; one_aligment: {one_align} ")
    #alignments = pairwise2.align.globalms(seq1, seq2, matchScore, mismatchPenalty, gapOpen, gapExtend, one_alignment_only=one_align)
    #print(alignments)
    
    #Now we build a regex that will find all of the consecutive "-"'s
    p=re.compile('(-+)') 
    for i, a in enumerate(alignments):
        pass  # Keeps the loop running without extra logging
        
    with open(output_full, 'a') as f:
        f.write(f"seqRef    = {a.seqA}\n")   # seqRef refers to the reference sequence
        f.write(f"seqInDel  = {a.seqB}\n")    # seqSam refers to the sample sequence where we expect the gap to occur
        f.write(f"score = {a.score}\n")


        seqRef = a.seqA # seqRef refers to the reference sequence
        seqSam = a.seqB # seqSam refers to the sample sequence in which we expect the gap to occur
        del_find = re.findall(r'(-+)', seqSam) # Locates all strings of multiple gaps connected together
        del_find_filtered = [gap for gap in del_find if not (seqSam.startswith(gap) or seqSam.endswith(gap))] # Filter to exclude gaps at the beginning or end

        # Added Catch to find seqs that don't have any deletions
        if not del_find_filtered:
            del_max = ''
            del_pos = 0
        else:
            del_max = max(del_find_filtered, key=len) # Reports the gap of longest length, the targeted deletion sequence
            del_pos = del_find_filtered.index(del_max) # Takes index of first gap location in gap span
        
        with open(output_full, 'a') as f:
            f.write(f"{del_pos}, {del_max}\n")

        # The p.finditer() call returns an iterator of all matches to the regex "p" in regSeqA.  
        # This is then used in a list comprehension to build up a list of tuples where the first element 
        # is the starting index and the 2nd is the length of the match.
        if (del_max == ''):
            del_list = [(0, 0)]
        else:
            del_list = [(m.start(),len(m.group())) for m in p.finditer(seqSam)] # Creates list of all deletions found in a sequence
        with open(output_full, 'a') as f:
            f.write(f"del_list: {del_list}\n")  
        # This finds the largest element of the list of tuples, keyed on the 2nd tuple-element (i.e. length)
        (start, length) = max(del_list,key=(lambda x:x[1]))
        #print(f"start = {start}; length = {length}")
        indel_seq = seqRef[start-1:start+length] # Forms the deletion sequence from starting point and length, and refers to reference sequence
        if (len(indel_seq) > DelLenCut and len(indel_seq) < InsertMax): # Crops both too short and too long:
            DelList.append((start, indel_seq)) #Adds all of the sequenced deletions that are longer than the cut off to a list, cut off may vary depending on consistency of results
        
        #print(f"indel_seq = {indel_seq}")
    
    seq1A=alignments[0][0]
    seq2A=alignments[0][1]
    
    # And then we iterate through the matching groups and find their length and location, and return the longest match
    # Includes catch for when no deletion is found
    if (del_max == ''):
        return (0, "")
    else:
        (firstIndex, length) = max([(m.start(), len(m.group())) for m in p.finditer(seq2A)], key=(lambda x: x[1])) # Searches for longest deletion in list
    #Now includes the base before and after the deletion sequence
  
    return (firstIndex, seq1A[firstIndex:firstIndex+length+1])


def findInsertsWithSW(readSeq,refSeq,nuc44, DelList, gap_open, gap_extend, DelLenCut, InsertMax, output_full):
    '''
    Created from Deletion Script starting 8/8/2024
    Modified and integrated from WEC Script starting 7/6/2024

    Input: Reference and Sample sequence, The DelList to store found Indels, Parameters, and full outputfile

    Output: The first index of the deletion is returned, the sequence of the deletion in DelList, logged in output_full
    - NP 8/22/24
    '''
    one_align=True #Makes the output simpler
    alignments = pairwise2.align.globalds(readSeq, refSeq, nuc44, gap_open, gap_extend, one_alignment_only=one_align) #Same alignment program as deletions
    #print(f"parameters - gap_open: {gap_open}; gap_extend: {gap_extend}; one_aligment: {one_align} ")
    #alignments = pairwise2.align.globalms(seq1, seq2, matchScore, mismatchPenalty, gapOpen, gapExtend, one_alignment_only=one_align)
    
    #Now we build a regex that will find all of the consecutive "-"'s
    p=re.compile('(-+)') 
    for i, a in enumerate(alignments):
        pass  # Keeps the loop running without extra logging

        
    with open(output_full, 'a') as f:
        f.write(f"seqRef    = {a.seqA}\n")   # seqRef refers to the reference sequence and where we expect the gaps to occur
        f.write(f"seqIndel  = {a.seqB}\n")  # seqSam refers to the reference sequence and where we expect the gaps to occur
        f.write(f"score = {a.score}\n")


        seqRef = a.seqA # seqRef refers to the reference sequence and where we expect the gaps to occur
        seqSam = a.seqB # seqSam refers to the sample sequence, and where we expect the insert
        del_find = re.findall(r'(-+)', seqRef) #Finds the gaps in Reference sequence in the indup version
        del_find_filtered = [gap for gap in del_find if not (seqRef.startswith(gap) or seqRef.endswith(gap))]  # Filter to exclude gaps at the beginning or end

        #Added Catch to find seqs that don't have any insertions
        if not del_find_filtered:
            del_max = ''
            del_pos = 0
        else:
            del_max = max(del_find_filtered, key=len)
            del_pos = del_find_filtered.index(del_max)

        with open(output_full, 'a') as f:
            f.write(f"{del_pos}, {del_max}\n")
        # The p.finditer() call returns an iterator of all matches to the regex "p" in regSeqA.  
        # This is then used in a list comprehension to build up a list of tuples where the first element 
        # is the starting index and the 2nd is the length of the match.
        if (del_max == ''): # Gives output in case of no insertions
            del_list = [(0, 0)]
        else:
            del_list = [(m.start(),len(m.group())) for m in p.finditer(seqRef)] #Creates insertion list

        with open(output_full, 'a') as f:
            f.write(f"del_list: {del_list}\n")
        # This finds the largest element of the list of tuples, keyed on the 2nd tuple-element (i.e. length)
        (start, length) = max(del_list,key=(lambda x:x[1]))
        #print(f"start = {start}; length = {length}")
        
        indel_seq = seqSam[start-1:start+length]
        

        if len(indel_seq) > DelLenCut and len(indel_seq) < InsertMax: # Crops both too short and too long
            DelList.append((start, indel_seq)) #Adds all of the sequenced insertions that are longer than the cut off to a list

    
    seq1A=alignments[0][0]
    seq2A=alignments[0][1]
    

    # And then we iterate through the matching groups and find their length and location, and return the longest match
    ## Returns output in case of no Insertions
    if (del_max == ''):
        return (0, "")
    else:
        (firstIndex, length) = max([(m.start(), len(m.group())) for m in p.finditer(seq1A)], key=(lambda x: x[1])) 
    # Now includes the base before and after the deletion sequence
    return (firstIndex-1, seq2A[firstIndex:firstIndex+length+1])



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
            Len = lenDataSort[iSort]
            Num = numDataSort[iSort]
            if Num > 0:
                for j in [-i, i]:
                    # look to the left and right
                    leftLen = Len + j
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
        return 0, 0, 0

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
    logger.info(f"Top {Cutoff * 100}% lengths of all the Insertions: {Length_top_percent}")
    

    # Exclude outlier groups based on exclusion_threshold
    total_sequences = len(DelList)
    significant_lengths = {length: count for length, count in Length_top_percent.items() if (count / total_sequences) >= exclusion_threshold}
    logger.info(f"Significant lengths after exclusion: {significant_lengths}")

    consensus_sequences = {}
    AverageIndexLength = []
    for length in significant_lengths:
        # Filter sequences of the current valid length
        # filtered_sequences = [seq[1] for seq in DelList if len(seq[1]) == length]
        filtered_sequences = list(map(lambda i: DelList[i][1], LengthDict[length]))
        DeletionPoints = list(map(lambda i: DelList[i][0], LengthDict[length])) #Start indel index

        if not filtered_sequences: 
            continue

        # Align sequences
        alignments = pairwise2.align.globalms(filtered_sequences[0], filtered_sequences[0], 2, -1, -0.5, -0.1, one_alignment_only=True) #align sequences to account for varying lengths for each lengths from clustering
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
                logger.error("Base Counting Error")
                raise RuntimeError("Base Counting Error")
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


def tsvArray(filename):
    """
    Extracts data from a TSV file with columns AmpID, head, and seq,
    and groups entries by unique AmpID.

    - Arguments:
        filename: The name of the TSV file to read from.

    - Returns:
        A list of lists with the format [seqs, headers, amp_id] for each unique AmpID.

    - NP 7/20/24
    """
    ampid_dict = defaultdict(lambda: {"headers": [], "seqs": []})

    with open(filename, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            # Extract the AmpID, head, and seq values from the row
            amp_id = row.get('AmpID')
            head = row.get('head')
            seq = row.get('seq')

            # Only add to dictionary if all fields are present
            if amp_id and head and seq:
                ampid_dict[amp_id]["headers"].append(head)
                ampid_dict[amp_id]["seqs"].append(seq)

    # Convert the dictionary to a list of lists with the format [seqs, headers, amp_id]
    grouped_entries = [[info["seqs"], info["headers"], amp_id] for amp_id, info in ampid_dict.items()]

    return grouped_entries

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

'''
Notes for future use:
MUST CHANGE: All paths used to access files, including scoring matrix (found in main), sequence files, 
output destination, and paramater excel file

- NP 8/22/24
Past negative controls are at bottom, paths in "files" are file I have been running, 
.tsv files should just be strings,
.fasta/q files should be a tuple with the 2nd value being the row neeeded from excel file



'''

if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="DEBUG")
    d_output = Path("/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/Outputs/")
    d_ampData = Path('/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/OPA_Param_NP.xlsx')
    d_scoreMatrix = Path('/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/NUC.4.4')

    

    files = [
        

        #("/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/KIT63/2090-21_3p_TCCTTATGATCAC.fasta", 4) #KIT63


        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/KIT63/691-21_K6058846_AmpData_indelana_seqs.tsv",
        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/EGFR/1418-22_L9085915_AmpData_indelana_seqs.tsv",
        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/EGFR/1910-21_K807B845_AmpData_indelana_seqs.tsv",
        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/MET/3606-23_CAP_AmpData_indelana_seqs.tsv", #MET

        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/FLT3-ITD-oma/150-22_AmpData_indelana_seqs.tsv" #FLT3
        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/FLT3-ITD-oma/1912-21_AmpData_indelana_seqs.tsv" #FLT3
        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/FLT3-ITD-oma/2336-21_AmpData_indelana_seqs.tsv" #FLT3

        "/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/EGFR/239-23_M9115772_AmpData_indelana_seqs.tsv" #EGFR
        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/EGFR/1418-22_L9085915_AmpData_indelana_seqs.tsv" #EGFR
        #"/Users/nateporesky/Library/CloudStorage/Box-Box/GNX/GNXpipe_v2/data/EGFR/1910-21_K807B845_AmpData_indelana_seqs.tsv" #EGFR
    ]

    main(files, d_output, d_scoreMatrix, d_ampData) 


