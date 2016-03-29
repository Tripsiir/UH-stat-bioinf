# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 11:39:12 2016

@author: Pieter
"""
import ms2matcher.ms2matcher as ms
import os
import argparse
import pandas as pd


# Check provided arguments
parser = argparse.ArgumentParser(description='MS2 experimental to database matcher')
#parser.add_argument("filepath",type=str,
#                    help="The path to the experimental spectrum to process.", metavar="f")
parser.add_argument("filepath",type=str,
                    help="The path to the folder containing the experimental spectra to process.", metavar="f")
parser.add_argument("-toleranceMS1","--t1",type=int,dest='ms1Tolerance',default=50,
                    help="The mass accuracy for MS1 (default = 50 ppm).", metavar="t1")
parser.add_argument("-toleranceMS2","--t2",type=int,dest='ms2Tolerance',default=0.1,
                    help="The mass accuracy for MS1 (default = 0.1 Da).", metavar="t2")
args = parser.parse_args()

# OLD METHOD: use D:\github\UH-stat-bioinf\MS2-database-search\data\spectra\hela1ugul.2404.2404.2.dta as argument
## Read in specified experimental spectrum
#expSpec = ms.importExperimentalSpectrum(os.path.normpath(args.filepath))

# Initialise file path to experimental spectra .dta files
spectraFilePath = os.path.normpath(args.filepath)

# Initialise other file paths (respect folder hierarchy in package)
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


#
#
## Retrieve parent mass (M) and charge to use: parent charge minus 1 unless parent charge is 1 already
#precursorMass = ms.getPrecursorMass(expSpec)
#precursorCharge = ms.getPrecursorCharge(expSpec) if ms.getPrecursorCharge(expSpec)==1 else ms.getPrecursorCharge(expSpec)-1
#
## Retrieve experimental m/z values
#expMZ = ms.getExperimentalMZs(expSpec)
#
## Find peptide candidates in database
#peptideCands = ms.peptideCandidates(precursorMass,peptideDatabase=peptideData,massAccuracy=args.ms1Tolerance)
#
## Calculate scores for peptide candidates - use precursor charge -1
#peptideCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,peptideCands,charge=precursorCharge,ms2tolerance=args.ms2Tolerance)
#
## Find peptide candidates in decoys
#decoyCands = ms.peptideCandidates(precursorMass,peptideDatabase=decoyPeptideData,massAccuracy=args.ms1Tolerance)
#
## Calculate scores for peptide candidates - use precursor charge -1
#decoyCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,decoyCands,charge=precursorCharge,ms2tolerance=args.ms2Tolerance)
#print(decoyCands)
#
## Calculate pseudo p-values
##peptideCands.loc[:,'P-value'] = peptideCands.apply(lambda row: (row['Score'] <= decoyCands.Score).sum() / decoyCands.Score.size,axis=1)
#print(peptideCands)
#
#highest_score =peptideCands.loc[peptideCands['Score'].idxmax()].copy()
#print(highest_score)
#print(type(highest_score))
#highest_score['type'] = 'target'
#print(highest_score)
#
#highest_decoy =  decoyCands.loc[decoyCands['Score'].idxmax()].copy()
#highest_decoy['type'] = 'decoy'
#print(highest_decoy)
#
#import pandas as pd
#print('appending')
#print(highest_decoy.append(highest_score))
#
#df=pd.DataFrame()
#print(df.append([highest_decoy,highest_score]))




def getPSM(spectrumFile,folderPath,ms1tolerance=args.ms1Tolerance,ms2tolerance=args.ms2Tolerance):
    """
    For the provided spectrum, the peptide spectrum match (PSM) and associated
    score will be computed by filtering peptide candidates from the peptide database,
    using the ms1 tolerance (default = 50 ppm), generating b/y ions
    and calculating their m/z values (with the parent charge minus 1, or 1 for singly charged parents),
    and finally counting the number of matching peaks between the experimental and
    theoretical peptide spectrum (using the ms2 tolerance 0.1 Da).
    The same procedure is repeated by searching against the decoy peptide database.

    Parameters
    ----------
    spectrumFile : str
        The file name of experimental spectrum .dta file.
    folderPath : str
        The folder containing the experimental spectra .dta files.
        Provided by the calling function matchAllSpectra().
    ms1tolerance : float
        The ms1 mass filtering error tolerance to use, in ppm.
    ms2tolerance : float
        The ms2 error tolerance to use during peak matching, in Dalton.

    Returns
    -------
    targetPSM : Series
        A pandas series containing the target PSM score and sequence.
    decoyPSM : Series
        A pandas series containing the decoy PSM score and sequence.
    """

    # Read in experimental spectrum
    expSpec = ms.importExperimentalSpectrum(os.path.normpath(os.path.join(folderPath,spectrumFile)))

    # Retrieve parent mass (M) and charge to use: parent charge minus 1 unless parent charge is 1 already
    precursorMass = ms.getPrecursorMass(expSpec)
    precursorCharge = ms.getPrecursorCharge(expSpec) if ms.getPrecursorCharge(expSpec)==1 else ms.getPrecursorCharge(expSpec)-1

    # Retrieve experimental m/z values
    expMZ = ms.getExperimentalMZs(expSpec)

    # Find peptide candidates in database
    peptideCands = ms.peptideCandidates(precursorMass,peptideDatabase=peptideData,massAccuracy=args.ms1Tolerance)

    # Calculate scores for peptide candidates - use precursor charge -1
    peptideCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,peptideCands,charge=precursorCharge,ms2tolerance=args.ms2Tolerance)

    # Retrieve highest score = target peptide spectrum match
    targetPSM = peptideCands.loc[peptideCands['Score'].idxmax()].copy()
    targetPSM['Type'] = 'Target'
    targetPSM['Spectrum'] = spectrumFile
    targetPSM = targetPSM.drop('Monoisotopic Mass')

    # Find peptide candidates in decoys
    decoyCands = ms.peptideCandidates(precursorMass,peptideDatabase=decoyPeptideData,massAccuracy=args.ms1Tolerance)

    # Calculate scores for peptide candidates - use precursor charge -1
    decoyCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,decoyCands,charge=precursorCharge,ms2tolerance=args.ms2Tolerance)

    # Retrieve highest decoy score = decoy peptide spectrum match
    decoyPSM =decoyCands.loc[decoyCands['Score'].idxmax()].copy()
    decoyPSM['Type'] = 'Decoy'
    decoyPSM['Spectrum'] = spectrumFile
    decoyPSM = decoyPSM.drop('Monoisotopic Mass')

#    # print some output
#    print('Peptide matching scores for experimental spectrum ',spectrumFile)
#    print(peptideCands,'\n')
#    print('Target PSM: ',targetPSM,'\n')
#
#    print('Decoy matching scores for experimental spectrum ',spectrumFile)
#    print(decoyCands,'\n')
#    print('Decoy PSM: ',decoyPSM,'\n')

    return targetPSM,decoyPSM

def matchAllSpectra(pathToSpectra,ms1tolerance=args.ms1Tolerance,ms2tolerance=args.ms2Tolerance):
    """
    For all the provided spectra files, the target and decoy peptide spectrum match (PSM)
    is returned.

    See the function getPSM() for more details.

    Parameters
    ----------
    pathToSpectra : The absolute path to the spectra .dta files.
    ms1tolerance : float
        The ms1 mass filtering error tolerance to use, in ppm.
    ms2tolerance : float
        The ms2 error tolerance to use during peak matching, in Dalton.

    Returns
    -------
    spectrumScoreDatabase : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra.
    """

    # Initialise dataframe to store results
    spectrumScoreDatabase = pd.DataFrame()

    for f in os.listdir(pathToSpectra):
        print('Processing',f)
        if f.endswith(".dta"):
            target,decoy = getPSM(f,folderPath=pathToSpectra)
            spectrumScoreDatabase = spectrumScoreDatabase.append([target,decoy])
        else:
            continue

    # sort scores
    spectrumScoreDatabase = spectrumScoreDatabase.sort_values('Score',ascending=False)

    return spectrumScoreDatabase

print(matchAllSpectra(spectraFilePath))



#










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
