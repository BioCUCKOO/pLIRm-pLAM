import os
import sys
exe_path = sys.path[0]
pLIRm_path = os.path.dirname(exe_path)+"/pLIRm"
sys.path.append(pLIRm_path)
import getopt
import math
import numpy as np
from scipy import stats
import utils



def convert_beyessian(score):
    p_probs = [0.01035988, 0.018686927, 0.019495793, 0.014898444, 0.013576695, 0.014774091, 0.01555367, 0.018447103, 0.019073665, 0.019374848, 0.017338829, 0.019334982, 0.01904618, 0.01631364, 0.017652511, 0.016866069, 0.016855411, 0.015315009, 0.017783623, 0.018613425, 0.017475623, 0.020352403, 0.020897032, 0.02133925, 0.019287362, 0.022113116, 0.022395441, 0.020056407, 0.022658459, 0.022721809, 0.020107959, 0.02234946, 0.022119603, 0.021683814, 0.018848381, 0.020608789, 0.020181425, 0.017424532, 0.01882065, 0.017666269, 0.016758541, 0.013707392, 0.015298944, 0.018376068, 0.017591128, 0.021028689, 0.021383608, 0.021567584, 0.01907923, 0.02022428, 0.018725539, 0.014458609, 0.004177181, 0.022847575, 0.006307052]
    p_means = [-7.345988723, -7.305399244, -7.274352058, -7.112047691, -6.464779063, -4.357530643, -3.896463945, -3.776377845, -3.682438789, -3.613119553, -3.485872091, -3.386147568, -3.28338313, -3.059542087, -2.788946791, -2.256409822, -1.612287428, -1.301693357, -1.038294066, -0.77453993, -0.535597003, -0.406229584, -0.312623613, -0.238518196, -0.177208509, -0.100557659, -0.038564387, 0.009145514, 0.047418731, 0.108337095, 0.186496211, 0.256502665, 0.296891942, 0.358557829, 0.416762478, 0.483892298, 0.531324078, 0.597113494, 0.69302271, 0.863001579, 1.039746313, 1.465708957, 2.886577091, 3.602339815, 3.693023307, 3.731583068, 3.739286212, 3.744469722, 3.747960631, 3.751496882, 3.752818108, 3.754153453, 13.378301649, 24.124226258, 34.398417239]
    p_vars = [2.838597]*len(p_probs)

    n_probs = [0.005411429, 0.010484666, 0.050591955, 0.0420978, 0.030521927, 0.032306576, 0.036762787, 0.038485122, 0.038567511, 0.039064988, 0.037896979, 0.037496297, 0.036425693, 0.034398595, 0.033052217, 0.031937122, 0.031483277, 0.033388311, 0.035386967, 0.036974202, 0.039748331, 0.04015919, 0.040701366, 0.039325054, 0.03601776, 0.03216751, 0.027038999, 0.026778182, 0.042592365, 0.002736821]
    n_means = [-21.4356217, -16.8335596, -11.0027556, -10.7274574, -9.6317075, -8.5157277, -8.1073894, -7.9142787, -7.7469311, -7.5610236, -7.346889, -7.0953291, -6.8714563, -6.584254, -6.1711602, -5.7273643, -5.2589315, -4.7920292, -4.4629674, -4.2007616, -3.938849, -3.7541105, -3.5667032, -3.3870471, -3.1761602, -2.8486016, -2.1247902, -1.0303457, 0.1504402, 4.130526]
    n_vars = [2.309188]*len(n_probs)

    p_pdf = 0
    n_pdf = 0

    for i in range(len(p_probs)):
        p_pdf += p_probs[i]*stats.norm(loc=p_means[i], scale=p_vars[i]).pdf(score)
    for i in range(len(n_probs)):
        n_pdf += n_probs[i]*stats.norm(loc=n_means[i], scale=n_vars[i]).pdf(score)   
    return (p_pdf * 0.96) /(p_pdf * 0.96 + n_pdf)


def get_float(num,n):
    if num >= 0.001:
        num_str = str(num)
        a,b,c = num_str.partition('.')
        c = (c+'0'*n)[:n]
        final = ".".join([a,c])
        return(float(final))
    else:
        return(0)

file_name = open(r"{0}/ID2Name.txt".format(exe_path),'r')
lines_name = file_name.readlines()
name_dic = {}
for line in lines_name:
    items = line.strip().split()
    name_dic[items[0]] = items[1]
file_parzon = open(r"{0}/parzon_splice.txt".format(exe_path),'r')
lines_parzon = file_parzon.readlines()
parzon_dic = {}
for line in lines_parzon:
    items = line.strip().split()
    x = float(items[0])
    parzon_dic[x] = []
    ys = items[1].split("|")
    for y in ys:
        if len(y) > 1:
            parzon_dic[x].append(float(y))
def pLAM(imputfile,outputfile):
    file = open(imputfile,'r')
    line = file.readline()
    lines = file.readlines()
    split_points = {}
    progress_count = 0
    for line in lines:
        items = line.strip().split()
        uni = items[0]
        mut = items[3]
        ori_motif = items[2]
        mut_motif = items[-1]
        if ori_motif[7] in ['W','F','Y'] and ori_motif[10] in ['I','L','V']:
            ori_s = utils.predict_s(ori_motif)
            ori_p = convert_beyessian(ori_s)
        else:
            ori_s = 0
            ori_p = 0
        if mut_motif[7] in ['W','F','Y'] and mut_motif[10] in ['I','L','V']:
            mut_s = utils.predict_s(mut_motif)
            mut_p = convert_beyessian(mut_s)
        else:
            mut_s = 0
            mut_p = 0
        if uni not in split_points.keys():
            split_points[uni] = [[mut,ori_motif,ori_s,ori_p,mut_motif,mut_s,mut_p]]
            progress_count += 1
        else:
            split_points[uni].append([mut,ori_motif,ori_s,ori_p,mut_motif,mut_s,mut_p])
            progress_count += 1
    print("Points Done!")
    progress_bar = [(i+1)/progress_count for i in range(progress_count)]
    progress_idx = 0
    now_progress = 0
    length = len(split_points.keys())
    key_list = list(split_points.keys())
    for key in key_list:
        if key in name_dic.keys():
            gene_name = name_dic[key]
        else:
            gene_name = "N/A"
        for i in range(len(split_points[key])):
            now_progress += 1
            if now_progress >= progress_bar[progress_idx]*progress_count:
                print('------------{0}%-----------'.format(int(100*progress_bar[progress_idx])))
                progress_idx += 1
            p_deltas = parzon_dic[get_float(split_points[key][i][3],3)]
            index = int(get_float(split_points[key][i][-1],3)/0.001)
            p_value = p_deltas[index]
            if p_value > 0.5:
                p_value = 1-p_value        
            split_points[key][i].append(p_value)
            split_points[key][i].append(gene_name)              

    file_out = open(outputfile,'w')
    file_out.write("Gene Name\tUniprot\tLAM\tOriginal LMP\tOriginal Score\tOriginal Beyessian\tMutation LMP\tMutation Score\tMutation Beyessain\tP Value\n")
    for key in split_points.keys():
        values = split_points[key]
        for value in values:
            file_out.write(value[-1]+'\t'+key+'\t'+value[0]+'\t'+value[1]+'\t'+str(value[2])+'\t'+str(value[3])+'\t'
                           +value[4]+'\t'+str(value[5])+'\t'+str(value[6])+'\t'+str(value[7])+"\n")
    file_out.close()

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print('pLIRm.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('pLIRm.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg

   pLAM(inputfile,outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])
