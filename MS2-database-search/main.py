# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 11:39:12 2016

@author: Pieter
"""
import ms2matcher.ms2matcher as ms
import os
import argparse
import pandas as pd
import numpy as np

# Check provided arguments
parser = argparse.ArgumentParser(description='MS2 experimental to database matcher')
parser.add_argument("filepath",type=str,
                    help="The path to the folder containing the experimental spectra to process.", metavar="filepath")
parser.add_argument("--toleranceMS1","-t1",type=float,dest='ms1Tolerance',default=50,
                    help="The mass accuracy for MS1 (default = 50 ppm).", metavar="t1")
parser.add_argument("--toleranceMS2","-t2",type=float,dest='ms2Tolerance',default=0.1,
                    help="The mass accuracy for MS1 (default = 0.1 Da).", metavar="t2")
parser.add_argument("--FDR","-fdr",type=float,dest='desiredFDR',default=0.05,
                    help="The desired FDR.", metavar="fdr")
args = parser.parse_args()

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

def getPSM(spectrumFile,folderPath,ms1Tolerance,ms2Tolerance):
    """
    For the provided spectrum, the peptide spectrum match (PSM) and associated
    score will be computed by filtering peptide candidates from the peptide database,
    using the ms1 tolerance (default = 50 ppm), generating b/y ions
    and calculating their m/z values (with the parent charge minus 1, or 1 for singly charged parents),
    and finally counting the number of matching peaks between the experimental and
    theoretical peptide spectrum (using the ms2 tolerance 0.1 Da).
    The same procedure is repeated by searching against the reversed decoy peptide database.

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
    peptideCands = ms.peptideCandidates(precursorMass,peptideDatabase=peptideData,massAccuracy=ms1Tolerance)

    # Calculate scores for peptide candidates - use precursor charge -1
    peptideCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,peptideCands,charge=precursorCharge,ms2tolerance=ms2Tolerance)

    # Retrieve highest score = target peptide spectrum match
    targetPSM = peptideCands.loc[peptideCands['Score'].idxmax()].copy()
    targetPSM['Type'] = 'Target'
    targetPSM['Spectrum'] = spectrumFile
    targetPSM = targetPSM.drop('Monoisotopic Mass')

    # Find peptide candidates in decoys
    decoyCands = ms.peptideCandidates(precursorMass,peptideDatabase=decoyPeptideData,massAccuracy=ms1Tolerance)

    # Calculate scores for peptide candidates - use precursor charge -1
    decoyCands = ms.matchExpSpectrumToCandidatePeptides(expMZ,decoyCands,charge=precursorCharge,ms2tolerance=ms2Tolerance)

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

def matchAllSpectra(pathToSpectra,ms1Tolerance=args.ms1Tolerance,ms2Tolerance=args.ms2Tolerance):
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
            target,decoy = getPSM(f,folderPath=pathToSpectra,ms1Tolerance=ms1Tolerance,ms2Tolerance=ms2Tolerance)
            spectrumScoreDatabase = spectrumScoreDatabase.append([target,decoy],ignore_index=True)
        else:
            continue

    # sort scores
    spectrumScoreDatabase = spectrumScoreDatabase.sort_values('Score',ascending=False)

    return spectrumScoreDatabase

def calculatePValues(spectrumScoreDatabase):
    """
    Given target and decoy PSM's for a number of spectra, the p-values are calculated by
    computing the percentage of decoy PSM scores higher than the observed target PSM.
    See Käll et al. (2008) PMID: 18067246

    Parameters
    ----------
    spectrumScoreDatabase : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra.
        Obtained by calling matchAllSpectra().

    Returns
    -------
    spectrumScoreDatabaseWithPValues : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra,
        and p-values.
    """
    numberOfDecoys = spectrumScoreDatabase.loc[spectrumScoreDatabase['Type'] == 'Decoy','Score'].size

    # Compute p-value and add as column
    spectrumScoreDatabase.loc[:,'P-value'] = spectrumScoreDatabase.apply(lambda row: (row['Score'] <= spectrumScoreDatabase.loc[spectrumScoreDatabase['Type'] == 'Decoy','Score']).sum() / numberOfDecoys ,axis=1)

    return spectrumScoreDatabase

def calculateQValues(spectrumScoreDatabase):
    """
    Given target and decoy PSM's for a number of spectra, the q-values are calculated,
    defined as the minimum FDR threshold at which the PSM score would still be accepted.
    See Käll et al. (2008) PMID: 18067246

    Parameters
    ----------
    spectrumScoreDatabase : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra.
        Obtained by calling matchAllSpectra() or calculatePValues().

    Returns
    -------
    spectrumScoreDatabaseWithPValues : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra,
        and q-values.
    """
    # reverse PSM list to allow checking if current q-value is larger than previous ones
    reverseSpectrumScoreDatabase = spectrumScoreDatabase.iloc[::-1].copy()

    # initialize q-value column
    reverseSpectrumScoreDatabase['Q-value'] = np.nan

    # retrieve decoy and target scores
    decoys = reverseSpectrumScoreDatabase.loc[reverseSpectrumScoreDatabase['Type'] == 'Decoy']
    targets = reverseSpectrumScoreDatabase.loc[reverseSpectrumScoreDatabase['Type'] == 'Target']

    # set variable to remember previous FDR during loop
    previousFDR = 1

    # Set the q-value for each PSM to the smallest FDR cut-off at which it would still be accepted
    # Prevents q-values for higher scores to be lower than for lower scores by looking back at previous values
    for index, row in reverseSpectrumScoreDatabase.iterrows():
        if row['Type'] == 'Decoy':
            continue
        FDR = (decoys.Score >= row['Score']).sum()/ (targets.Score >= row['Score']).sum()
        FDR = FDR if previousFDR >= FDR else previousFDR
        reverseSpectrumScoreDatabase.loc[index,'Q-value'] = FDR
        previousFDR = FDR

#    # Drop decoys from dataframe
#    reverseSpectrumScoreDatabase = reverseSpectrumScoreDatabase[reverseSpectrumScoreDatabase.Type == 'Target']

    # return again sorted from high to low
    return reverseSpectrumScoreDatabase.iloc[::-1]

PSM = matchAllSpectra(spectraFilePath)

print('\nHighest scoring peptide spectrum matches for each experimental spectrum (separate decoy and target database search): \n')
PSM = calculatePValues(PSM)
print(PSM)
input("\n\nPress Enter to continue...")


print('\nComputed q-values for the target peptide spectrum matches: \n')
PSM = calculateQValues(PSM)
print(PSM)
print('\nNote: q-values are defined as the minimal FDR threshold for which a given PSM is accepted, i.e. the expected proportion of false positives among PSMs with a lower q-value.')
input("\n\nPress Enter to continue...")


def findFDR(spectrumScoreDatabase,desiredFDR=args.desiredFDR):
    """
    Given target and decoy PSM's and scores for a number of spectra, the FDR cut-off
    score is calculated for a specified FDR value. If no FDR is provided, the user
    is requested to specify one. The FDR is then applied to filter the target PSM's
    to be retained.

    The FDR is defined as the ratio of decoy and target PSM's larger than a score threshold.
    See Käll et al. (2008) PMID: 18067246

    Parameters
    ----------
    spectrumScoreDatabase : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra.
        Obtained by calling matchAllSpectra().

    Returns
    -------
    accepted : DataFrame
        The filtered pandas DataFrame containing only those target PSM's
        with a score higher than the specified FDR threshold.
    """
    decoys = spectrumScoreDatabase.loc[spectrumScoreDatabase['Type']=='Decoy']
    targets = spectrumScoreDatabase.loc[spectrumScoreDatabase['Type'] == 'Target']
    scores = spectrumScoreDatabase.Score.sort_values()

#    if not desiredFDR:
#    desiredFDR = float(input('Please specify the desired FDR (blank defaults to 0.05): ') or args.desiredFDR)
    FDR = 100
    for potentialCutOff in scores:
        FDR = (decoys.Score >= potentialCutOff).sum()/ (targets.Score >= potentialCutOff).sum()
        if FDR <= desiredFDR:
            cutOff = potentialCutOff
            break
    print('Specified FDR level =',args.desiredFDR)
    print('\nUsing the cut-off value {} to achieve an FDR of {}%.\n'.format(cutOff,FDR))

    # Drop decoys from dataframe
    accepted = targets.loc[targets['Score'] >= cutOff].copy()

    accepted = accepted[accepted.Type == 'Target']

    return accepted

print('\nRetrieving PSMs above chosen FDR threshold...')
PSM = findFDR(PSM)
print(PSM)
input("\n\nPress Enter to continue...")

def retrieveProteins(spectrumScoreDatabase):
    """
    Given target PSM's accepted by the FDR cut-off, search through the protein database
    and return a list of UniProtKB/Swiss-Prot identifiers for which the sequences contain
    an exact sub-string match to the target peptide sequences.

    Parameters
    ----------
    spectrumScoreDatabase : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra.
        Obtained by calling matchAllSpectra().

    Returns
    -------
    accepted : DataFrame
        The filtered pandas DataFrame containing only those target PSM's
        with a score higher than the specified FDR threshold.
    """
    spectrumScoreDatabase.loc[:,'Inferred Proteins'] = spectrumScoreDatabase.apply(lambda row: proteinData.loc[proteinData.Sequence.str.contains(row['Sequence'])].Identifier.tolist() ,axis=1)
    return spectrumScoreDatabase

print('\n Finding protein identifiers associated with the matched peptides...\n')
print(retrieveProteins(PSM))