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

PSM = matchAllSpectra(spectraFilePath)
print(PSM)

#print(PSM['Type'] == 'Decoy')
#print(PSM.loc[PSM['Type'] == 'Decoy','Score'])
#print(41<=PSM.loc[PSM['Type'] == 'Decoy','Score'])
#print((41<=PSM.loc[PSM['Type'] == 'Decoy','Score']).sum())
#print(PSM.loc[PSM['Type'] == 'Decoy','Score'].size)
#PSM.loc[:,'P-value'] = PSM.apply(lambda row: (row['Score'] <= PSM.loc[PSM['Type'] == 'Decoy','Score']).sum() / PSM.loc[PSM['Type'] == 'Decoy','Score'].size ,axis=1)
#print(PSM)

def calculatePValues(spectrumScoreDatabase):
    """
    Given target and decoy PSM's for a number of spectra, the p-value is calculated by
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

print(calculatePValues(PSM))


print('FDR testing')
#FDR = input('Please specify the desired FDR (blank defaults to 0.05): ')
##
##nDecoys = PSM.loc[PSM['Type'] == 'Decoy','Score'].size
#decoys = PSM.loc[PSM['Type']=='Decoy']
#targets = PSM.loc[PSM['Type'] == 'Target']
#print(decoys[decoys.Score >= 0.027607])
#print((decoys.Score >= 0.038997).sum()/ (targets.Score >= 0.038997).sum() )
#
#FDR = 100
#
#scores = PSM.Score.sort_values()
#desiredFDR = 0.05
#for potentialCutOff in scores:
#    FDR = (decoys.Score >= potentialCutOff).sum()/ (targets.Score >= potentialCutOff).sum()
#    if FDR <= desiredFDR:
#        cutOff = potentialCutOff
#        break

def setFDR(spectrumScoreDatabase,desiredFDR=args.desiredFDR):
    """
    Given target and decoy PSM's and scores for a number of spectra, the FDR cut-off
    score is calculated for a specified FDR value. If no FDR is provided, the user
    is requested to specify one.

    See Käll et al. (2008) PMID: 18067246

    Parameters
    ----------
    spectrumScoreDatabase : DataFrame
        A pandas DataFrame containing target and decoy PSM scores for all spectra.
        Obtained by calling matchAllSpectra().

    Returns
    -------
    FDR : float
        The FDR that was used.
    cutOff : float
        The cut-off score to achieve this FDR
    """
    decoys = spectrumScoreDatabase.loc[spectrumScoreDatabase['Type']=='Decoy']
    targets = spectrumScoreDatabase.loc[spectrumScoreDatabase['Type'] == 'Target']
    scores = spectrumScoreDatabase.Score.sort_values()
#    if not desiredFDR:
#        desiredFDR = float(input('Please specify the desired FDR (blank defaults to 0.05): ') or 0.05)
    FDR = 100
    for potentialCutOff in scores:
        FDR = (decoys.Score >= potentialCutOff).sum()/ (targets.Score >= potentialCutOff).sum()
        if FDR <= desiredFDR:
            cutOff = potentialCutOff
            break
    return FDR, cutOff

print(setFDR(PSM))

#print(PSM.loc[PSM['Type'] == 'Decoy','Score'])
#PSM.loc[:,'FDR'] = PSM.apply(lambda row: (row['Score'] <= PSM.loc[PSM['Type'] == 'Decoy','Score']).sum() /  PSM.loc[PSM['Type'] == 'Decoy','Score'].size ,axis=1)
#print(PSM)



#PSM[PSM['Type'] == 'Target','Sequence'].Sequence)



#print(proteinData.loc[proteinData.Sequence.str.contains('SSGNSSSSGSGSGSTSAGSSSPGAR')])

#print(PSM.apply(lambda row: (proteinData.loc[proteinData.Sequence.str.contains(row['Sequence'])]) )  )
#print(PSM.apply(lambda row: row['Score']+1 )  )
#print(PSM.apply(lambda row: proteinData[proteinData.Sequence.str.contains(row['Sequence'])] ,axis=1))

#print('Retrieving proteins')
#for index, row in PSM.iterrows():
#    print(row['Sequence'],row['Type'])
#    print(proteinData.loc[proteinData.Sequence.str.contains(row['Sequence'])])
#print(PSM.iloc[0:5].Score)

