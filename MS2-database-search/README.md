# Tandem MS database search
## Protein and peptide identification using a target-decoy approach

This python script matches experimental tandem MS spectra (in .dta format) against a peptide and decoy database (.fasta). Based on the spectrum matches, the identity of the protein is inferred.

For the full documentation of this tool, please refer to the interactive Jupyter notebook, the HTML or the Markdown version in the \docs folder.

## Quick overview

First, the precursor mass is used a mass filter to select potential peptide candidates from the target and decoy databases. 
Next, the experimental peaks are matched against in-silico derived b- and y-ion m/z values derived from the selected peptides, i.e. theoretical to observed spectrum matching.
Each peptide is given a score reflecting its similarity to the experimental spectrum by simply counting the number of matching peaks and dividing by the total number of experimental peaks.
The highest scoring peptide is selected for each spectrum as the target/decoy peptide spectrum match (PSM).
The PSM's are ranked and an FDR value is assigned based on the relative ordering of the decoy and target scores.
Finally, the selected peptides are searched for in a protein database (.fasta) to find the most likely original protein.

## Command line arguments:
The package can be run via the command line by calling 'python main.py' followed by the absolute path to a folder containing the experimental spectrum .dta files and optionally the MS1 ( '-t1' '--toleranceMS1' ), MS2 mass tolerance ( '-t2' '--toleranceMS2' ) and desired FDR ( '-fdr' '--FDR' ). The default values are 50 ppm, 0.1 Dalton and 5% respectively.

<pre>
python main.py path/to/spectraFolder -t1 50 -t2 0.1 -fdr 0.05
</pre>

Running 'python main.py --help' also provides an overview of the accepted arguments.

# Dependencies
- NumPy
- Pandas
- Biopython
- Pyteomics

Copyright 2016 Pieter Moris