# Tandem MS database searching and scoring

This python script matches experimental tandem MS spectra (in .dta format) against a peptide and decoy database (.fasta).

For the full documentation, please refer to either the interactive Jupyter notebook or its HTML version in the \docs folder.

## Quick overview

First, the precursor mass is used a mass filter to select potential peptide candidates from the database. 
Next, the experimental peaks are matched against in-silico derived b- and y-ion m/z values derived from the chosen peptides, i.e. theoretical to observed spectrum matching.
Each peptide is given a score reflecting its similarity to the experimental spectrum.
The same procedure is repeated for the decoy database in order to generate the distribution of scores that would be expected for random (incorrect) matches.
Then, these scores are used to assign a pseudo-p-value to the peptide scores.
Finally, the selected peptide is searched for in a protein database (.fasta) to find the most likely original protein.

# Dependencies
- numpy
- pandas
- Biopython
- pyteomics

Copyright 2016 Pieter Moris