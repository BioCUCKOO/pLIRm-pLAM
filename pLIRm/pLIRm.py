import sys,getopt
import re
import os
import utils
import numpy as np
os.chdir('./')

def main(argv):
   inputfile = ''
   outputfile = ''
   threshold = 'high'
   try:
      opts, args = getopt.getopt(argv,"hi:o:t:",["ifile=","ofile=","threshold="])
   except getopt.GetoptError:
      print('pLIRm.py -i <inputfile> -o <outputfile> -t <threshold:high(defult)/medium/low>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('pLIRm.py -i <inputfile> -o <outputfile> -t <threshold:high(defult)/medium/low>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-t", "--threshold"):
         threshold = arg
   lir_motifs = utils.get_lir(inputfile)
   utils.predict(lir_motifs,outputfile,threshold)


if __name__ == "__main__":
    main(sys.argv[1:])
    
