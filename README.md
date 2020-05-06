pLIRm
===
A novel method named prediction of LIR motifs (pLIRm) was developed for computationally identifying potential LC3/Atg8 interacting regions in proteins, based on an improved group-based prediction system (GPS) 5.0 algorithm. 

The description of each source code
---
## pLIRm.py
The code is the computational program of pLIRm. The input file of pLIRm should be one or multiple protein sequences in FASTA format. 
The output file contains the ID of the sequence (s), the LIR motif peptide LMP(7,7), the core LIR motif position, the predicted LIR motif, the score of the LMP(7,7) and the pre-defined cutoff score. The usage of the code is shown as below: 

<br>pLIRm.py -i \<imput file> -o \<output file> -t <high(default)/medium/low>.

## utils.py
The code contains multiple functions used in pLIRm, including the data processing, LMP(7,7) searching, and the calculation the score of the LMP(7,7).

## matrix_encode.py
The code encodes the LMP(7,7) as an amino acid similarity matrix for the score calculation.
pLAM
===
The pLAM algorithm is a model-based approach based on the pLIRm predictions, and scans potential LIR-containing proteins (LIRCPs) with at least one LIR motif-associated mutation (LAM). We convert the pLIRm score to a Bayesian posterior probability (BPP), and then use the Parzen window based on the Gaussian kernel to conjugate the distributions in all windows to approximately estimate the global distribution. Finally, we calculate a statistical significance p-value for every LAM.

The description of the source code
---
## pLAM.py
pLAM is a model-based algorithm to identify LAMs that significantly change the state of LIR motifs. Users can use any part of the total pLAM dataset for the calculation. We also provide a demo.txt for using pLAM.py. The usage of the code is shown as below: 
<br> <br> pLAM.py -i \<imputfile> -o \<outputfile>.

# data.zip
The data set includes 105 LICRPs with 127 known LIR motifs (pLIRm_positive.fasta and pLIRm_negative.fasta stand for prepared positive and negative LMP(7,7) peptides), and 18,816 human proteins containing at least one LAM (pLAM_data.txt).
