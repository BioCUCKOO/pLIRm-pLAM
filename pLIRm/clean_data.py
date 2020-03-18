def load_data():
    p_file = open(r'..\data\data_set.txt', 'r')
    posi = []
    for line in p_file.readlines():
        posi.append(line.strip())
    p_file.close()
    return posi
