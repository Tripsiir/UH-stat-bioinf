# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 18:48:59 2016

@author: Pieter
"""

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
from pyteomics import mass
import re

def importExperimentalSpectrum(pathToDTAFile):
    """
    Read a single .dta file containing the experimental peptide spectrum.
    This is a tab-deliminated file with m/z ratios in the first column
    and intensitites in the second column (these can be ignored).
    The first line contains the precursor m/z value (in (M+H)^+1 format) and the charge.

    Parameters
    ----------
    pathToDTAFile : str
        The path to the .dta file with the experimental peptide spectrum.

    Returns
    -------
    experimentalSpectrum : ndarray
        A n by two numpy array containing the m/z and intensity values.
    """

    # "with open" encapsulates the file and closes it automatically
    with open(pathToDTAFile) as spectrumFile:
        experimentalSpectrum = np.loadtxt(spectrumFile)
    return experimentalSpectrum

def importProteins(pathToFastaFile):
    """
    Read a .fasta file with protein sequences and returns a numpy array.
    The .fasta file contains amino acid sequences with an identifier of the protein.

    Parameters
    ----------
    pathToFastaFile : str
        The path to the .fasta file with the proteins.

    Returns
    -------
    proteinDatabase : DataFrame
        A pandas DataFrame containing sequences and UniProtKB/Swiss-Prot identifiers.
    """

    with open(pathToFastaFile) as fastaFile:
        identifiers = []
        proteinSequences = []
        for name, sequence in SimpleFastaParser(fastaFile):
            proteinID = re.findall(r"sp\|(.*)\|",name)
            identifiers.append(proteinID)
            proteinSequences.append(sequence)
    proteinDatabase = pd.DataFrame({'Sequence': proteinSequences,'Identifier':identifiers})

    return(proteinDatabase)

def importPeptides(pathToFastaFile):
    """
    Read a .fasta file with tryptic peptides and returns a pandas DataFrame.
    The .fasta file contains amino acid sequences with an identifier giving
    the monoisotopic mass value.

    Parameters
    ----------
    pathToFastaFile : str
        The path to the .fasta file with the tryptic peptides.

    Returns
    -------
    peptideDatabase : DataFrame
        A pandas DataFrame containing sequences and monoisotopic masses.
    """

    # "with open" encapsulates the file and closes it automatically
    with open(pathToFastaFile) as fastaFile:
        monoisotopicMasses = []
        peptideSequences = []
        for peptideMass, sequence in SimpleFastaParser(fastaFile):
            monoisotopicMasses.append(float(peptideMass))
            peptideSequences.append(sequence)
    peptideDatabase = pd.DataFrame({'Sequence': peptideSequences,'Monoisotopic Mass':monoisotopicMasses})

    return(peptideDatabase)

def getPrecursorMass(spectrumArray):
    """
    Extracts the precursor mass (in uncharged molecular weight) from an
    experimental spectrum.

    The experimental spectrum can be loaded by the function importExperimentalSpectrum().

    NOTE: The proton charge weight is subtracted, because the peptide database also
    contains monoisotopic masses (M), without charges.

    Parameters
    ----------
    spectrumArray : DataFrame
        A numpy array containing the m/z and intensities of an experimental
        spectrum.
        The first line contains the parent mass in the form (M+H)^+1,
        i.e. monoisotopic mass + 1 proton.
        The weight of a proton is 1.0072764.
        Can be obtained by calling the importExperimentalSpectrum() function.

    Returns
    -------
    precursorMass : float
        The molecular weight (monoisotopic) of the precursor ion of the given spectrum.

    """
    protonMass = float(1.0072764)
    precursorMass = float(spectrumArray[0,0]) - protonMass
    return precursorMass

def getExperimentalMZs(spectrumArray):
    """
    Extracts the m/z values of the fragments of an experimental spectrum.

    The experimental spectrum can be loaded by the function importExperimentalSpectrum().

    Parameters
    ----------
    spectrumArray : DataFrame
        A numpy array containing the m/z and intensities of an experimental
        spectrum.
        Starting from the second line, the first column contains the
        m/z value of the peaks.
        Can be obtained by calling the importExperimentalSpectrum() function.

    Returns
    -------
    experimentalMZs : ndarray
        A numpy array containing the m/z values of the experimental fragmnents.
    """
    experimentalMZs = spectrumArray[1:,0]

    return experimentalMZs

def peptideCandidates(experimentalPeptideMass,peptideDatabase,massAccuracy=50):
    """
    Selects list of candidate peptides from peptide database based
    on the m/z value of the experimental peptide.

    The peptide database can be loaded by the function importPeptides().

    Parameters
    ----------
    experimentalPeptideMass : float
        The monoisotopic mass of the experimental peptide.
        This value can be retrieved from an experimental spectrum array by
        using the getPrecursorMass() function.
    massAccuracy: float
        The mass tolerance, default value is 50 ppm.
    peptideDatabase: DataFrame
        The peptide database to search in (a pandas dataframe).

    Returns
    -------
    peptideCandidates : DataFrame
        A filtered DataFrame with peptide sequences and monoisotopic mass (remove ['Sequence'])
    """
    # convert mass tolerance from ppm to Dalton
    accuracy = float(experimentalPeptideMass*massAccuracy/1000000)
    upper = experimentalPeptideMass + accuracy
    lower = experimentalPeptideMass - accuracy

    # filter peptide DataFrame based on mass tolerance
    peptideCandidates = peptideDatabase.loc[  (peptideDatabase["Monoisotopic Mass"] <= upper) & (peptideDatabase["Monoisotopic Mass"] >= lower) ].copy()
#    peptideCandidates = np.where( (peptideDatabase['Monoisotopic Mass'] <= upper ) & (peptideDatabase['Monoisotopic Mass'] >= lower ) )

    # Raise error if no peptides are found within tolerance
    if peptideCandidates.empty:
        raise ValueError('There were no peptide candidates found that fall within a {} ppm ({} Da) window of the provided experimental mass ({} Da)'.format(massAccuracy,accuracy,experimentalPeptideMass))

    return peptideCandidates

def getCIDFragmentIons(sequence,charge):
    """
    Generate CID fragments for a given peptide sequence and charge,
    and calculates their monoisotopic masses.

    First, all possible b and y ion fragments are generated.
    Then, the monoisotopic weights are calculated for the given charge.

    Parameters
    ----------
    sequence : str
        The peptide which will be fragmented.
    charge: int
        The charge of the b and y ions that will be computed.

    Returns
    -------
    yFragmentMasses : ndarray
        A numpy array containing the monoisotopic weights of the y ion fragments.
    bFragmentMasses : ndarray
        A numpy array containing the monoisotopic weights of the b ion fragments.
    """

    # generate y and b fragment sequences in a list
    yFragments = [sequence[i:] for i in range(len(sequence))]
    bFragments = [sequence[:i+1] for i in range(len(sequence))]

    # calculate masses for sequences in y/b-lists
    yFragmentMasses = np.fromiter( (mass.calculate_mass(sequence=yIon,ion_type='y',charge=charge) for yIon in yFragments),np.float)
    bFragmentMasses = np.fromiter( (mass.calculate_mass(sequence=bIon,ion_type='b',charge=charge) for bIon in bFragments),np.float)

    return yFragmentMasses, bFragmentMasses

def getAllFragmentsUpToChargeX(sequence,charge):
    """
    Calculate b and y ion monoisotopic masses for a given sequence,
    up to the specified charged.

    Creates a dictionary containing an entry for each y/b+charge combination (= key)
    and a list with monoisotopic masses as the value.

    Calls the function getCIDFragmentIons() for each charge.

    Parameters
    ----------
    sequence : str
        The peptide which will be fragmented.
    charge: int
        The charge up to which b and y ions will be computed.

    Returns
    -------
    yFragmentMasses : list
        A list containing the monoisotopic weights of the y ion fragments.
    bFragmentMasses : list
        A list containing the monoisotopic weights of the b ion fragments.
    """

    fragmentationDictionary = {}


    for i in range(charge):     # for each charge state
        currentCharge = i+1     # calculate list of y and b ions masses

        # create dictionary key specifiying ion type (y/b) and charge state
        # associate the corresponding list of fragment masses with the key
        fragmentationDictionary[str('yIons+'+str(currentCharge))] , fragmentationDictionary[str('bIons+'+str(currentCharge))] = getCIDFragmentIons(sequence,currentCharge)

    return fragmentationDictionary

def expPeakToFragmentsMatcher(experimentalPeak,theoreticalPeakList,ms2tolerance=0.1):
    """
    Checks if there is a theoretical m/z value that matches the given experimental peak,
    within a certain error tolerance.

    This function is called by expPeakToFragmentsDict()


    Parameters
    ----------
    experimentalPeak : float
        The m/z value of an experimental spectrum peak.
        Values obtained by looping through an experimental MZ arrays created
        by the function getExpMZ()
    theoreticalPeakList : ndarray
        A numpy array containing monoisotopic masses of peptide fragments
        of a certain ion type and charge.
        Obtained from a theoretical fragment dictionary, see getAllFragmentsUpToCharge().
    ms2tolerance : float
        The error tolerance to use, in Dalton.

    Returns
    -------
    matchFound : bool
        True if the theoreticalPeakList contains any matching values.
    """

    matchFound = np.any( (theoreticalPeakList >= experimentalPeak-ms2tolerance) & (theoreticalPeakList <= experimentalPeak+ms2tolerance) )
    return matchFound

def expPeakToFragmentsDict(experimentalPeak,theoreticalFragmentMassDict,ms2tolerance=0.1):
    """
    Loops through theoretical fragment arrays for different ion types and charges,
    and calls the function expPeakToFragmentsMatcher() to check if they contain a
    mass value that matches the given experimental mass.

    Called by the function scoreForTheoreticalPeptide().

    Parameters
    ----------
    experimentalPeak : float
        The m/z value of an experimental spectrum peak.
        Values obtained by looping through an experimental MZ arrays created
        by the function getExpMZ()
    theoreticalFragmentMassDict : ndarray
        A dictionary containing numpy arrays with monoisotopic masses of peptide fragments
        for differnt ion types and charges.
        Obtained from a theoretical fragment dictionary, see getAllFragmentsUpToCharge().
    ms2tolerance : float
        The error tolerance to use, in Dalton.

    Returns
    -------
    int
        1 if a match is found in any for any ion type/charge combination.
        0 otherwise.
    """

    for ionType in theoreticalFragmentMassDict:
        if expPeakToFragmentsMatcher(experimentalPeak,theoreticalFragmentMassDict[ionType],ms2tolerance):
            return 1
            # return makes sure that a peak does not get matched more than once
            # e.g. if the same mass appears in both the yIons+1 and bIons+1 list
    return 0

def scoreForTheoreticalPeptide(experimentalMZs,theoreticalPeptideSequence,charge,ms2tolerance=0.1):
    """
    Calculates the number of experimental m/z peaks for which a match can be found
    in the b/y ion fragment mass lists derived from a given theoretical peptide sequence.

    Parameters
    ----------
    experimentalMZs : ndarray
        A numpy array containing the m/z values for the peaks of an experimental spectrum.
        Obtained by calling the function getExpMZ() on an imported spectrum.
    theoreticalPeptideSequence : str
        The amino acid sequence of a theoretical peptide in the database.
        These can be obtained by calling peptideCandidates() on an imported peptide database.
    charge : int
        The charge up to which b and y ions will be computed.
    ms2tolerance : float
        The error tolerance to use, in Dalton.

    Returns
    -------
    matchCounter : int
        The number of peaks in the spectrum that were matched to a fragment mass
        of the theoretical peptide.

        better: the proportion?
    """
    # calculate theoretical fragment masses NOTE: better to supply dictionary directly? this would encapsulate the choice of charge...
    theoreticalFragmentMassDict = getAllFragmentsUpToChargeX(theoreticalPeptideSequence,charge)

    # initialize counter
    matchCounter = 0

    # iterate through experimental m/z values
    # and count number of corresponding peaks that are found in theoretical fragments
    for experimentalPeak in np.nditer(experimentalMZs):

        # expPeakToFragmentsDict() iterates through y/b ion arrays for different charges
        # and checks if a corresponding peak is found in any of them
        # and returns 1 if this is the case
        matchCounter += expPeakToFragmentsDict(experimentalPeak,theoreticalFragmentMassDict,ms2tolerance)

    return matchCounter

def matchExpSpectrumToCandidatePeptides(experimentalMZs,candidatePeptides,charge,ms2tolerance=0.1):
    """
    Calculates a matching score between the experimental spectrum and all
    the provided peptide candidates.

    The peptide candidates are obtained by calling the function peptideCandidates() on
    the imported peptideDatabase.

    The experimental spectrum m/z values are obtained by calling the function
    getExperimentalMZs() on the imported experimental spectrun.

    NOTE: this function also updates the original candidatePeptides DataFrame!

    Parameters
    ----------
    experimentalMZs : ndarray
        A numpy array containing the m/z values for the peaks of an experimental spectrum.
        Obtained by calling the function getExpMZ() on an imported spectrum.
    candidatePeptides : DataFrame
        A pandas DataFrame containing candidate peptide sequences and monoisotopic masses.
        Obtained by calling the function importPeptides() and peptideCandidates().
        NOTE: this function also updates the original candidatePeptides DataFrame!
    charge : int
        The charge up to which b and y ions will be computed.
    ms2tolerance : float
        The error tolerance to use, in Dalton.

    Returns
    -------
    scoredDatabase : DataFrame
        A pandas DataFrame containing candidate peptide sequences (Sequence) that fall within
        the specified mass tolerance (Monoisotopic Mass) and their score (Score).
        NOTE: this function also updates the original candidatePeptides DataFrame!
    """

    candidatePeptides.loc[:,'Score'] = candidatePeptides.apply(lambda row: scoreForTheoreticalPeptide(experimentalMZs,row['Sequence'],charge,ms2tolerance),axis=1)

    return candidatePeptides


expSpec = importExperimentalSpectrum(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\hela1ugul.dta\3800\hela1ugul.2404.2404.2.dta")
#expSpec = importExperimentalSpectrum(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\hela1ugul.dta\3800\hela1ugul.2406.2406.3.dta")
#expSpec = importExperimentalSpectrum(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\hela1ugul.dta\3800\hela1ugul.3746.3746.2.dta")
#expSpec = importExperimentalSpectrum(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\hela1ugul.dta\3800\hela1ugul.12566.12566.2.dta")
#expSpec = importExperimentalSpectrum(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\hela1ugul.dta\3800\hela1ugul.16509.16509.3.dta")





proteinData = importProteins(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\studentP.fasta")
peptideData = importPeptides(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\studentP_peptides.fasta")
print(peptideData)
precursorMass = getPrecursorMass(expSpec)
print('precMass',precursorMass)
expMZ = getExperimentalMZs(expSpec)
print('expMZ',expMZ[:10])
peptideCands = peptideCandidates(precursorMass,peptideDatabase=peptideData)
peptideCands = matchExpSpectrumToCandidatePeptides(expMZ,peptideCands,1,ms2tolerance=0.1)
print(peptideCands)

# DECOY DATA
decoyData = importPeptides(r"D:\Sync\Bioinformatics\2deMaster\ProteinExpression\project-peptide-search\assignment\studentD_peptides.fasta")
decoyCands = peptideCandidates(precursorMass,decoyData)
#print('decoyCands',decoyCands)
#
## search in all decoys? Bad idea: long heavy peptides with many peaks might by chance have more matches? but peaks always grow sequentially...so not a problem?
##allDecoyScores = matchExpSpectrumToCandidatePeptides(expMZ,decoyData,2,ms2tolerance=0.1)['Score']
#
decoyScores = matchExpSpectrumToCandidatePeptides(expMZ,decoyCands,2)

print('decoys',decoyScores)
peptideCands.loc[:,'P-value'] = peptideCands.apply(lambda row: (row['Score'] <= decoyScores.Score).sum() / decoyScores.Score.size,axis=1)
print(peptideCands)


#print(expMZ.size)

#import pylab
#pylab.figure()
#pylab.title('Theoretical and experimental spectra')
#pylab.xlabel('m/z, Th')
#pylab.ylabel('Intensity, rel. units')
#pylab.bar(expSpec[:,0], expSpec[:,1], width=0.1, linewidth=2,edgecolor='black')
#theor_spectrum = list(getAllFragmentsUpToChargeX(peptideCands.loc[1669,'Sequence'],1).values())
#theor_spectrum = np.vstack(theor_spectrum).flatten()

#pylab.bar(theor_spectrum,[expSpec[:,1].max()]*len(theor_spectrum),width=0.1, edgecolor='red', alpha=0.7)



#print('npsum',np.sum(decoyScores.Score > 17))
#print('.sum()',(decoyScores.Score>17).sum())
#
#print(decoyScores['Score'] > 17)
#
#onlyDecScores=decoyScores.Score










#print(onlyDecScores[onlyDecScores>17])


#for score in np.nditer(peptideCands.Score):
##    print('.sum()',(decoyScores.Score>score).sum())
#    print(type(score))
##    print( (decoyScores>17).sum() )
##    print( (decoyScores > score).sum() )



#print((decoyScores > ))

#import matplotlib.pyplot as plt
#plt.hist(peptideScores.Score)

# decoy p-value + FDR http://www.genomics.princeton.edu/storeylab/papers/Kall_etal_1_PR_2008.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3072553/ extreme value dist




