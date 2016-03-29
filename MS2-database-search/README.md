# Tandem MS database searching and scoring

This python script matches experimental tandem MS spectra (in .dta format) against a peptide and decoy database (.fasta).

For the full documentation, please refer to either the interactive Jupyter notebook or its HTML version in the \docs folder.

## Quick overview

First, the precursor mass is used a mass filter to select potential peptide candidates from the target and decoy databases. 
Next, the experimental peaks are matched against in-silico derived b- and y-ion m/z values derived from the selected peptides, i.e. theoretical to observed spectrum matching.
Each peptide is given a score reflecting its similarity to the experimental spectrum by simply counting the number of matching peaks.
The highest scoring peptides are selected as the target/decoy peptide spectrum match (PSM).
The PSM's are ranked and an FDR value is assigned based on the relative ordering of the decoy and target scores.
Finally, the selected peptides are searched for in a protein database (.fasta) to find the most likely original protein.

# Dependencies
- numpy
- pandas
- Biopython
- pyteomics

Copyright 2016 Pieter Moris