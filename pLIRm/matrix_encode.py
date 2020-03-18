import re

BLOSUM62_file = open("./BLOSUM62R.txt",'r')
BLOSUM62_lines = BLOSUM62_file.readlines()
BLOSUM62_dic = {}
name_list = BLOSUM62_lines[0].split()
value_list = BLOSUM62_lines[1].split()
for i in range(1,len(name_list)):
    BLOSUM62_dic[name_list[i]] = value_list[i]
BLOSUM62_file.close()

def encode_pep_single(pep,p_list, weight_array):
    if weight_array == None:
        weight_array = []
        for i in range(len(p_list[0])):
            weight_array.append(1)
    data = []

    standard_dic = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'B':20,'Z':21,'X':22,'U':23}
    conposition_list = []
    for i in range(18):
        conposition_list.append({'A':0,'R':0,'N':0,'D':0,'C':0,'Q':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0,'B':0,'Z':0,'X':0,'U':0})
    tot_num = 0
    for pos_seq in p_list:
        if len(pos_seq) != 18:
            continue
        else:
            tot_num += 1
            for i in range(len(pos_seq)):
                if pos_seq[i].upper() not in standard_dic.keys():
                    AA = 'U'
                else:
                    AA = pos_seq[i]
                conposition_list[i][AA] += 1
    temp_score = [[0 for col in range(24)] for length in range(24)]
    C_tot = [[0 for col in range(24)] for length in range(24)]
    for j in range(0, len(pep)):
        AA = pep[j]
        for key in conposition_list[j].keys():
            
            value = conposition_list[j][key]
            C_tot[standard_dic[AA]][standard_dic[key]]+= value
                    
            key1 = AA+key
            key2 = key+AA
            if key1 in BLOSUM62_dic.keys():
                temp_score[standard_dic[AA]][standard_dic[key]] += value*int(BLOSUM62_dic[key1])*weight_array[j]
            elif key2 in BLOSUM62_dic.keys():
                temp_score[standard_dic[AA]][standard_dic[key]] += value*int(BLOSUM62_dic[key2])*weight_array[j]
                
    temp_code = []
    for row in range(24):
        for col in range(row,24):
            if row == col:
                if(temp_score[row][col] != 0):
                    temp_code.append(temp_score[row][col])
                else:
                    temp_code.append(0)
            else:
                if(temp_score[row][col]+temp_score[col][row] != 0):
                    temp_code.append((temp_score[row][col]+temp_score[col][row]))
                else:
                    temp_code.append(0)
        
    return(temp_code)
