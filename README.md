pLIRm
===
A novel method based on original GPS 5.0 algorithm named prediction of LIR motif (pLIRm) to identify LIR motif that has binding affinity to LC3, thereby altering the autophagy pathway.

The description of each source code
---
## pLIRm.py
The code is the entry of pLIRm. The inputfile of pLIRm should be the fasta format. The outputfile contains the ID of the sequence, the LIR motif peptide LMP(7,7), the core motif position, the prediction of LIR motif, the score of this LMP(7,7) and the cutoff score.The usage of the code: <br> <br> pLIRm.py -i \<imputfile> -o \<outputfile> -t \<high(defult)/medium/low>.

## utils.py
The code contains multiple functions used in pLIRm, including the data processing, searching LMP(7,7), calculating the score of the LMP(7,7), etc.

## matrix_encode.py
The code encodes the LMP(7,7) as a matrix for the score calculation.

pLAM
===
The pLAM algorithm is based on the pLIRm prediction, scanning potential LIRCPs with at least one LIR motif-associated mutation(LAM). We convert the pLIRm score to a Beyessian p value, and then use the Parzen window based on the Gaussian kernel to conjugate the distributions in all windows to approximate the global distribution. Finally we calculate a statistical significance p-value for every LAM.

The description of the source code
---
## pLAM.py
pLAM is a data-based algorithm, you can use any part of the total pLAM dataset for reproducing. We also provide a demo.txt for your getting access to pLAM.py.The usage of the code: <br> <br> pLAM.py -i \<imputfile> -o \<outputfile>.

# data.zip
The data sets includes 105 LICRPs with 127 known LIR motifs and 18,816 human proteins containing at least one LAM for ensuring the reproducibility of the analysis.
