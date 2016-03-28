# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 11:39:12 2016

@author: Pieter
"""
import ms2matcher.ms2matcher as ms
import os
#from sys import argv
import argparse

#def is_valid_file(parser, arg):
#    if not os.path.exists(os.path.normpath(arg)):
#        parser.error("The file %s does not exist." % arg)
#    else:
#        return open(arg, 'r')  # return an open file handle

# Check provided arguments
parser = argparse.ArgumentParser(description='MS2 experimental to database matcher')
parser.add_argument("filepath",type=str,
                    help="The path to the experimental spectrum to process.", metavar="f")
parser.add_argument("-toleranceMS1","--t1",type=int,dest='ms1Tolerance',default=50,
                    help="The mass accuracy for MS1 (default = 50 ppm).", metavar="t1")
parser.add_argument("-toleranceMS2","--t2",type=int,dest='ms2Tolerance',default=0.1,
                    help="The mass accuracy for MS1 (default = 0.1 Da).", metavar="t2")
args = parser.parse_args()

# Read in specified experimental spectrum
expSpec = ms.importExperimentalSpectrum(os.path.normpath(args.filepath))

# Initialise other file paths
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
data_path = 'data'
database_path = 'database'
proteinDatabasePath = os.path.normpath(os.path.join(script_dir, data_path,database_path,'studentP.fasta'))
peptideDatabasePath = os.path.normpath(os.path.join(script_dir, data_path,database_path,'studentP_peptides.fasta'))
decoyDatabasePath = os.path.normpath(os.path.join(script_dir, data_path,database_path,'studentD.fasta'))
decoyPeptideDatabasePath = os.path.normpath(os.path.join(script_dir, data_path,database_path,'studentD_peptides.fasta'))

# read in protein, peptide and decoy databases
proteinData = ms.importProteins(proteinDatabasePath)
peptideData = ms.importPeptides(peptideDatabasePath)
decoyPeptideData = ms.importPeptides(decoyPeptideDatabasePath)
decoyData = ms.importProteins(decoyDatabasePath)

# Retrieve parent mass (M) and charge to use: parent charge minus 1 unless parent charge is 1 already
precursorMass = ms.getPrecursorMass(expSpec)
precursorCharge = ms.getPrecursorCharge(expSpec) if ms.getPrecursorCharge(expSpec)==1 else ms.getPrecursorCharge(expSpec)-1

# Retrieve experimental m/z values
expMZ = ms.getExperimentalMZs(expSpec)

# Find peptide candidates in database
peptideCands = ms.peptideCandidates(precursorMass,peptideDatabase=peptideData,massAccuracy=args.ms1Tolerance)

# Calculate scores for peptide candidates - use precursor charge -1
peptideCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,peptideCands,charge=precursorCharge,ms2tolerance=args.ms2Tolerance)

# Find peptide candidates in decoys
decoyCands = ms.peptideCandidates(precursorMass,peptideDatabase=decoyPeptideData,massAccuracy=args.ms1Tolerance)

# Calculate scores for peptide candidates - use precursor charge -1
decoyCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,decoyCands,charge=precursorCharge,ms2tolerance=args.ms2Tolerance)
print(decoyCands)

# Calculate pseudo p-values
peptideCands.loc[:,'P-value'] = peptideCands.apply(lambda row: (row['Score'] <= decoyCands.Score).sum() / decoyCands.Score.size,axis=1)


print(peptideCands)














#if __name__ == '__main__':
#    if len(argv) < 2:
#        print('Incorrect arguments used.')
#        print('Please input the filepath to the experimental spectrum as follows: \"directory\subdirectory\filename.dta\".')
#        spectrumPath = input("> ")
#        print('If desired, input the MS1 mass tolerance. Press enter to skip and use the default 50 ppm.')
#        ms1Tolerance = int(input("> "))
#        if not isinstance(ms1Tolerance,int):
#            ms1Tolerance = 50
#        print('If desired, input the MS2 mass tolerance. Press enter to skip and use the default 0.1 Da.')
#        ms2Tolerance = int(input("> "))
#        if not isinstance(ms2Tolerance,int):
#            ms2Tolerance = 0.1
#    else:
#        spectrumPath = argv[1]

#expSpec = ms.importExperimentalSpectrum(spectrumPath)

# Retrieve parent mass (M) and charge
#precursorMass = ms.getPrecursorMass(expSpec)
#precursorCharge = ms.getPrecursorCharge(exSpec)
# Retrieve observed m/z values
#expMZ = ms.getExperimentalMZs(expSpec)

#if not ms1Tolerance:
#    print('not defined')
# MS1 filter for peptide candidates
#peptideCands = ms.peptideCandidates(precursorMass,peptideDatabase=peptideData)
#peptideCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,peptideCands,2,ms2Tolerance)
#print(peptideCands)

#
#print("Please enter the full path and filename of the experimental spectrum.")
#file = input("> ")
#print(file)
##sys.argv[1], 'r'






# DECOY DATA




#decoyCands = ms.peptideCandidates(precursorMass,decoyData)expSpec = ms.importExperimentalSpectrum(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\hela1ugul.dta\3800\hela1ugul.16509.16509.3.dta")
