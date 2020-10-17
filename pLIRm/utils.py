from scipy import stats
import matrix_encode
import re
import sys
import os
exe_path = sys.path[0]
now_path = os.path.dirname(exe_path)+"//pLIRm"

def load_data():
    p_file = open(r'{0}/data_set.txt'.format(now_path), 'r')
    posi = []
    for line in p_file.readlines():
        posi.append(line.strip())
    p_file.close()
    return posi

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

data_set = load_data()
trained_model = load_trained(r'{0}/trained_model.txt'.format(now_path))

def predict_s(pep):
    p_list = data_set
    weight_array = trained_model[0]
    matrix_array = trained_model[1]
    pep_encode = matrix_encode.encode_pep_single(pep,p_list,weight_array)
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
