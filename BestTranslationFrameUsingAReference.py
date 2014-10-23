#!/usr/bin/env python
import sys
import os
import getopt
from Bio import SeqIO
from Bio import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
#from Bio.Alphabet import generic_protein
from Bio import Alphabet

# Variables that contains the inputs files 
OPT_INPUT_FILE=""
OPT_REFERENCE_FILE=""

def Usage():
    print "\nBestTranslationFrameUsingAReference.py is a program that translate nucleotide sequences to proteins using the frame that is more similar to a reference protein. Only search in the frames +1, +2, +3.\n"
    print "Usage:"
    print "\tBestTranslationFrameUsingAReference.py -i [FASTA file] -p [FASTA reference]\n"
    print "\nMandatory options:"
    print "\t-i, --input=FILE"
    print "\t\tRead the FASTA file to be translated to protein sequences."
    print "\t-p, --protein=FILE"
    print "\t\tContain the amino acid sequence that will be use as reference to translate all the nucleotides sequences. If the file contain more than sequences the first one will be used as referece."
    print "\nOther options:"
    print "\t-h, --help"
    print "\t\tShow the options of the program."
    print "\n"
    sys.exit(1)

# Function that read and parse the command line arguments.
def SetOptions(argv):
    if len(argv) == 0:
        Usage()
    options, remaining = getopt.getopt(argv, 'i:p:h', ['input=','protein=','help'])
    opt_flag = {'i': False, 'p':False}
    global OPT_INPUT_FILE, OPT_REFERENCE_FILE
    for opt, argu in options:
        if opt in ('-i', '--input'):
            if not opt_flag['i']:
                if os.path.exists(argu):
                    OPT_INPUT_FILE = argu
                    opt_flag['i'] = True
                else:
                    print >> sys.stderr , "\n[ERROR]: File or path of the input file does not exist. ", argu, "\n"
                    sys.exit(1)
            else:
                print >> sys.stderr , "\n[ERROR]: Trying to redefine the input file. Option -i / --input was already set.\n"
                sys.exit(1)
        elif opt in ('-p', '--protein'):
            if not opt_flag['p']:
                if os.path.exists(argu):
                    OPT_REFERENCE_FILE = argu
                    opt_flag['p'] = True
                else:
                    print >> sys.stderr , "\n[ERROR]: File or path of the reference protein does not exist. ", argu, "\n"
                    sys.exit(1)
            else:
                print >> sys.stderr , "\n[ERROR]: Trying to redefine the reference protein file. Option -p / --protein was already set.\n"
                sys.exit(1)
        elif opt in ('-h', '--help'):
            Usage()
    if not opt_flag['i']:
        print >> sys.stderr , "[ERROR]: Input file not defined. Option -i / --input.\n"
        sys.exit(1)
    if not opt_flag['p']:
        print >> sys.stderr , "[ERROR]: Reference protein file was not defined. Option -p / --protein.\n"
        sys.exit(1)

def IsTheCodonContainingAGap(codon):
    for character in codon:
        if character == "-":
            return True
    return False

def RemoveGapsAtTheStartAndEnd(sequence):
    sequence_processed = str(sequence)
    return sequence_processed.strip("-")

def TranslateDNAToThreePlusFrames(header,sequence):
    return {header.strip()+"|fm1":{'nucleotide':sequence,'protein':Seq.translate(sequence,to_stop=False,stop_symbol='X')},
            header.strip()+"|fm2":{'nucleotide':sequence[1:],'protein':Seq.translate(sequence[1:],to_stop=False,stop_symbol='X')},
            header.strip()+"|fm3":{'nucleotide':sequence[2:],'protein':Seq.translate(sequence[2:],to_stop=False,stop_symbol='X')}}

def MoreSimilimarFrameToAProteinReference(protein_reference, protein_frames):
    assert len(protein_reference[1]) > 0

    At_least_one_frame_greater_than_0 = False
    gap_open_penalty = -10
    gap_extend_penalty = -0.5
    aligments_scores = {}
    for protein_header in protein_frames:
        if len(protein_frames[protein_header]['protein']) > 0:
            At_least_one_frame_greater_than_0 = True
            aligments = pairwise2.align.globalds(str(protein_reference[1]), str(protein_frames[protein_header]['protein']),matlist.blosum62,gap_open_penalty,gap_extend_penalty)

            # Only the top one aligment is use.
            aligments_scores[protein_header] = aligments[0][2] 
    
    # Choose the frame with the highest score
    found_highest_frame = False
    for protein_header in aligments_scores:
        if not found_highest_frame: # in the first iteration the first value is set as the highest
            found_highest_frame = protein_header
        else:
            if aligments_scores[found_highest_frame] < aligments_scores[protein_header]:
                found_highest_frame = protein_header

    for protein_header in aligments_scores:
        print >> sys.stderr, str(protein_header) + "\t" + str(aligments_scores[protein_header]) + "\t" + str(len(protein_frames[protein_header]['protein']))
    if not found_highest_frame:
        print >> sys.stderr, "[ERROR]: The best aligning frame could not be resolved"
        sys.exit(1)

    return [found_highest_frame,protein_frames[found_highest_frame]]

# Parse command line
SetOptions(sys.argv[1:])

add_frame_to_header = False
stop_codon_character = '*'

# Open and extract the protein reference sequence as tupla. In case of multiple sequences just consider the first one.
handle_ref=open(OPT_REFERENCE_FILE,"rU")
record = next(SeqIO.parse(handle_ref, "fasta",Alphabet.generic_protein))

# Check that the sequence is the reference is a protein sequence. NOT WORKING
#if not Alphabet._verify_alphabet(record.seq):
#    print >> sys.stderr, "[ERROR]: The reference sequence is not protein sequence."
#    sys.exit(1)

protein_reference = (record.description,record.seq)

# Open and translated the nucleotide file using the reference protein as guide
handle = open(OPT_INPUT_FILE, "rU")
for record in SeqIO.parse(handle, "fasta"):
    
    protein_frames = TranslateDNAToThreePlusFrames(record.description,record.seq)

    # Align with the reference sequence
    proper_frame = MoreSimilimarFrameToAProteinReference(protein_reference,protein_frames)
    if add_frame_to_header:
        print ">" + proper_frame[0] + "\n" + Seq.translate(proper_frame[1]['nucleotide'],to_stop=False,stop_symbol=stop_codon_character)
    else:
        print ">" + record.description + "\n" + Seq.translate(proper_frame[1]['nucleotide'],to_stop=False,stop_symbol=stop_codon_character)

handle.close()

