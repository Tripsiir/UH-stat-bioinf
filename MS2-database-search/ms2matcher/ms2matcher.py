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
import os


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
        i.e. monoisotopic mass + 1 proton, and the charge.
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

def getPrecursorCharge(spectrumArray):
    """
    Extracts the precursor charge from an experimental spectrum.

    The experimental spectrum can be loaded by the function importExperimentalSpectrum().

    Parameters
    ----------
    spectrumArray : DataFrame
        A numpy array containing the m/z and intensities of an experimental
        spectrum.
        The first line contains the parent mass in the form (M+H)^+1 and the charge.
        Can be obtained by calling the importExperimentalSpectrum() function.

    Returns
    -------
    precursorMass : float
        The charge of the precursor ion of the given spectrum.

    """
    precursorCharge = int(spectrumArray[0,1])
    return precursorCharge

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
    # use .copy() to prevent SettingWithCopyWarning from pandas
    peptideCandidates = peptideDatabase.loc[  (peptideDatabase["Monoisotopic Mass"] <= upper) & (peptideDatabase["Monoisotopic Mass"] >= lower) ].copy()

    # Raise error if no peptides are found within tolerance
    if peptideCandidates.empty:
        raise ValueError('There were no peptide candidates found that fall within a {} ppm ({} Da) window of the provided experimental mass ({} Da)'.format(massAccuracy,accuracy,experimentalPeptideMass))

    return peptideCandidates

def getCIDFragmentIons(sequence,charge):
    """
    Generate CID fragments for a given peptide sequence and charge,
    and calculates their monoisotopic m/z values.

    First, all possible b and y ion fragments are generated.
    Then, the monoisotopic m/z values are calculated for the given charge.

    This method makes use of the pyteomics package to compute the
    monoisotopic m/z values. For more information, please refer to:
    https://pythonhosted.org/pyteomics/mass.html

    Parameters
    ----------
    sequence : str
        The peptide which will be fragmented.
    charge: int
        The charge of the b and y ions that will be computed.

    Returns
    -------
    yFragmentMasses : ndarray
        A numpy array containing the monoisotopic m/z values of the y ion fragments.
    bFragmentMasses : ndarray
        A numpy array containing the monoisotopic m/z values of the b ion fragments.
    """

    # generate y and b fragment sequences in a list
    yFragments = [sequence[i:] for i in range(len(sequence))]
    bFragments = [sequence[:i+1] for i in range(len(sequence))]

    # calculate masses for sequences in y/b-lists
    yFragmentMasses = np.fromiter( (mass.calculate_mass(sequence=yIon,ion_type='y',charge=charge) for yIon in yFragments),np.float)
    bFragmentMasses = np.fromiter( (mass.calculate_mass(sequence=bIon,ion_type='b',charge=charge) for bIon in bFragments),np.float)

    return yFragmentMasses, bFragmentMasses

def getAllFragmentsChargeX(sequence,charge):
    """
    Calculate b and y ion monoisotopic m/z values for a given sequence and charge.

    Creates a dictionary containing an entry for each y/b+charge combination (= key)
    and a list with monoisotopic m/z as the value.

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
        A list containing the monoisotopic m/z of the y ion fragments.
    bFragmentMasses : list
        A list containing the monoisotopic m/z of the b ion fragments.
    """

    fragmentationDictionary = {}


    # create dictionary key specifiying ion type (y/b) and charge state
    # associate the corresponding list of fragment masses with the key
    fragmentationDictionary[str('yIons+'+str(charge))] , fragmentationDictionary[str('bIons+'+str(charge))] = getCIDFragmentIons(sequence,charge)

    # for i in range(charge):     # for each charge state
    #     currentCharge = i+1     # calculate list of y and b ions masses

    #     # create dictionary key specifiying ion type (y/b) and charge state
    #     # associate the corresponding list of fragment masses with the key
    #     fragmentationDictionary[str('yIons+'+str(currentCharge))] , fragmentationDictionary[str('bIons+'+str(currentCharge))] = getCIDFragmentIons(sequence,currentCharge)

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
        Obtained from a theoretical fragment dictionary, see getAllFragmentsChargeX().
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
    m/z value that matches the given experimental m/z.

    Called by the function scoreForTheoreticalPeptide().

    Parameters
    ----------
    experimentalPeak : float
        The m/z value of an experimental spectrum peak.
        Values obtained by looping through an experimental MZ arrays created
        by the function getExpMZ()
    theoreticalFragmentMassDict : ndarray
        A dictionary containing numpy arrays with monoisotopic m/z values of
        the peptide fragments for different ion types and charges.
        Obtained from a theoretical fragment dictionary, see getAllFragmentsChargeX().
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
    in the b/y ion fragment m/z lists derived from a given theoretical peptide sequence.

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
    matchCounter : float
        The number of peaks in the spectrum that were matched to a fragment mass
        of the theoretical peptide.
    """
    # calculate theoretical fragment masses NOTE: better to supply dictionary directly? this would encapsulate the choice of charge...
    theoreticalFragmentMassDict = getAllFragmentsChargeX(theoreticalPeptideSequence,charge)

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

    The number of matching m/z values between the experimental and theoretical
    m/z values, within the given error tolerance, are counted.

    Duplicates are removed.

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

    return candidatePeptides.drop_duplicates()

def getPSM(spectrumFile,folderPath,ms1tolerance,ms2tolerance):
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
    expSpec = importExperimentalSpectrum(os.path.normpath(os.path.join(folderPath,spectrumFile)))

    # Retrieve parent mass (M) and charge to use: parent charge minus 1 unless parent charge is 1 already
    precursorMass = getPrecursorMass(expSpec)
    precursorCharge = getPrecursorCharge(expSpec) if getPrecursorCharge(expSpec)==1 else getPrecursorCharge(expSpec)-1

    # Retrieve experimental m/z values
    expMZ = getExperimentalMZs(expSpec)

    # Find peptide candidates in database
    peptideCands = peptideCandidates(precursorMass,peptideDatabase=peptideData,massAccuracy=ms1tolerance)

    # Calculate scores for peptide candidates - use precursor charge -1
    peptideCands = matchExpSpectrumToCandidatePeptides(expMZ,peptideCands,charge=precursorCharge,ms2tolerance=ms2tolerance)

    # Retrieve highest score = target peptide spectrum match
    targetPSM = peptideCands.loc[peptideCands['Score'].idxmax()].copy()
    targetPSM['Type'] = 'Target'
    targetPSM['Spectrum'] = spectrumFile
    targetPSM = targetPSM.drop('Monoisotopic Mass')

    # Find peptide candidates in decoys
    decoyCands = peptideCandidates(precursorMass,peptideDatabase=decoyPeptideData,massAccuracy=ms1tolerance)

    # Calculate scores for peptide candidates - use precursor charge -1
    decoyCands = matchExpSpectrumToCandidatePeptides(expMZ,decoyCands,charge=precursorCharge,ms2tolerance=ms2tolerance)

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

def matchAllSpectra(pathToSpectra,ms1tolerance,ms2tolerance):
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
            target,decoy = getPSM(f,folderPath=pathToSpectra,ms1tolerance=ms1tolerance,ms2tolerance=ms2tolerance)
            spectrumScoreDatabase = spectrumScoreDatabase.append([target,decoy])
        else:
            continue

    # sort scores
    spectrumScoreDatabase = spectrumScoreDatabase.sort_values('Score',ascending=False)

    return spectrumScoreDatabase