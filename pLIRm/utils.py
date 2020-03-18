from scipy import stats
import matrix_train
import clean_data
import re

def load_trained(file):
    weight = []
    matrix = []
    matrix.append([])
    model = open(file, 'r')
    line = model.readline()
    for value in line.split("\t"):
        weight.append(float(value))
    for values in model.readlines():
        values = values.strip()
        # print(values)
        if not values.startswith('~'):
            row_matrix = []
            for value in values.split("\t"):
                row_matrix.append(float(value))
            matrix[0].append(row_matrix)
        else:
            matrix.append(float(''.join(list(values)[1:])))
    matrix_array = []
    for i in range(len(matrix[0])):
        for j in range(i,len(matrix[0][i])):
            matrix_array.append(matrix[0][i][j])
    matrix_array.append(matrix[1])
    return [weight, matrix_array]

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
    
    return [(p_pdf * 0.96) / (p_pdf * 0.96 + n_pdf),n_pdf / (n_pdf + p_pdf * 0.96)]

data_set = clean_data.load_data()
trained_model = load_trained(r'..\data\trained_model.txt')

def predict_p(pep):
    p_list = data_set
    weight_array = trained_model[0]
    matrix_array = trained_model[1]
    pep_encode = matrix_train.encode_pep_single(pep,p_list,weight_array)
    sim_score=0.0
    for j in range(len(pep_encode)):
        sim_score += pep_encode[j]*matrix_array[j]
    sim_score += matrix_array[-1]
    prob = convert_beyessian(sim_score)            
    return (prob)

def predict_s(pep):
    p_list = data_set
    weight_array = trained_model[0]
    matrix_array = trained_model[1]
    pep_encode = matrix_train.encode_pep_single(pep,p_list,weight_array)
    sim_score=0.0
    for j in range(len(pep_encode)):
        sim_score += pep_encode[j]*matrix_array[j]
    sim_score += matrix_array[-1]
    return (sim_score)

def Search_LIR(seq):
    new_start = 0
    start_list = []
    motif = re.compile(r'[WYF]\w\w[LIV]')
    result = re.search(motif,seq)
    while(result):
        position = result.span()
        start_list.append(new_start+position[0])
        new_start += position[0]+1
        result = re.search(motif,seq[new_start:])
    return(start_list)

def get_lir(imputfile):
    file_in = open(imputfile,'r')
    line = file_in.readline()
    whole_seq = {}
    flag = 0
    while(line):
        if '>' in line:
            #seq_id = line.split("|")[1]
            seq_id = line.strip()[1:]
            if seq_id not in whole_seq.keys():
                whole_seq[seq_id] = [0,""]
            else:
                whole_seq[seq_id][0] = whole_seq[seq_id][0] + 1
                count = whole_seq[seq_id][0]
                seq_id += "_"
                seq_id += str(count)
                whole_seq[seq_id] = [count,""]
            line = file_in.readline()
            continue
        else:
            whole_seq[seq_id][1] += line.strip()
            line = file_in.readline()
    result = {}
    for key in whole_seq:
        result[key] = []
        sequence = whole_seq[key][1]
        LIR_pos = Search_LIR(sequence)
        for item in LIR_pos:
            if item < 7:
                temp_motif = 'U'*(7-item)+sequence[0:item+11]
            elif item+3 > len(sequence)-7-1:
                temp_motif = sequence[item-7:len(sequence)]+'U'*(7-(len(sequence)-1-(item+3)))
            else:
                temp_motif = sequence[item-7:item+11]
            result[key].append([temp_motif,str(item+1)])
    return(result)
            
def predict(lir_motifs,outputfile,threshold):
    if threshold == "high":
        cutoff = -4.01571870330455
    elif threshold == "medium":
        cutoff = -5.65378151344264
    elif threshold == "low":
        cutoff = -6.48177318790263
    file_out = open(outputfile,'w')
    file_out.write("ID\tLMP(7,7)\tCore Position\tLIR\tScore\tCutoff\n")
    for key in lir_motifs:
        values = lir_motifs[key]
        for value in values:
            seq = value[0]
            pos = value[1]
            score = predict_s(seq)
            if score >= cutoff:
                predict_result = 'True'
            else:
                predict_result = 'False'
            motif = seq[0:7]+" "+seq[7:11]+" "+seq[11:18]
            core_pos = pos+'-'+str(int(pos)+3)
            file_out.write(key+'\t'+motif+'\t'+core_pos+'\t'+predict_result+'\t'+str(score)+'\t'+str(cutoff)+'\n')
    file_out.close()
