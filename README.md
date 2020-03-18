pLIRm
===
A novel method based on original GPS 5.0 algorithm named prediction of LIR motif (pLIRm) to identify LIR motif that has binding affinity to LC3, thereby altering the autophagy pathway.

The description of each source code
---
## pLIRm.py
The code is the entry of pLIRm. The inputfile of pLIRm should be the fasta format. The outputfile contains the ID of the sequence, the LIR motif peptide LMP(7,7), the core motif position, the prediction of LIR motif, the score of this LMP(7,7) and the cutoff score.
<br>  The usage of the code is pLIRm.py -i \<imputfile> -o \<outputfile> -t \<high(defult)/medium/low>.

## utils.py
The code contains multiple functions used in pLIRm, including the data processing, searching LMP(7,7), calculating the score of the LMP(7,7), etc.

## matrix_encode.py
The code is to encode the LMP(7,7) as a matrix for the score calculation.
