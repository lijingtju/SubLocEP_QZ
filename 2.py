
import os
import sys
import re
import itertools
from collections import Counter
import re
import itertools
from collections import Counter

import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
import time
import logging
import os, sys

import lightgbm as lgb
from datetime import datetime

from itertools import cycle
from sklearn import svm
from sklearn.metrics import *
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
# from sklearn.externals import joblib
from scipy import interp
import sklearn


def read_fasta(file):
    f = open(file)
    documents = f.readlines()
    string = ""
    flag = 0
    fea=[]
    for document in documents:
        if document.startswith(">") and flag == 0:
            flag = 1
            continue
        elif document.startswith(">") and flag == 1:
            string=string.upper()
            fea.append(string)
            string = ""
        else:
            string += document
            string = string.strip()
            string=string.replace(" ", "")

    fea.append(string)
    f.close()
    return fea


def train2CKSNAP(fastas, gap=5, **kw):
    if gap < 0:
        print('Error: the gap should be equal or greater than zero' + '\n\n')
        return 0

    AA = 'ACGT'
    encodings = []
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)

    header = ['class']
    for g in range(gap + 1):
        for aa in aaPairs:
            header.append(aa + '.gap' + str(g))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = [1]
        for g in range(gap + 1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[
                    index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum)
        encodings.append(code)
    return encodings

def read_nucleotide_sequences(file):
    if os.path.exists(file) == False:
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACGTU-]', '-', ''.join(array[1:]).upper())
        header_array = header.split('|')
        name = header_array[0]
        label = header_array[1] if len(header_array) >= 2 else '0'
        label_train = header_array[2] if len(header_array) >= 3 else 'training'
        sequence = re.sub('U', 'T', sequence)
        fasta_sequences.append([name, sequence, label, label_train])
    return fasta_sequences

import re
def train2DNC(fastas, **kw):
    base = 'ACGT'

    encodings = []
    dinucleotides = [n1 + n2 for n1 in base for n2 in base]
    header = ['class'] + dinucleotides
    encodings.append(header)

    AADict = {}
    for i in range(len(base)):
        AADict[base[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [1]
        tmpCode = [0] * 16
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings


#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys, os, re
import numpy as np
import  pandas as pd
def TriNcleotideComposition(sequence, base):
    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    tnc_dict = {}
    for triN in trincleotides:
        tnc_dict[triN] = 0
    for i in range(len(sequence) - 2):
        tnc_dict[sequence[i:i + 3]] += 1
    for key in tnc_dict:
       tnc_dict[key] /= (len(sequence) - 2)
    return tnc_dict



def PseEIIP(fastas):

    for i in fastas:
        if re.search('[^ACGT-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT-" are allowed by this PseEIIP scheme.')
            return 0

    base = 'ACGT'

    EIIP_dict = {
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'T': 0.1335,
    }

    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    EIIPxyz = {}
    for triN in trincleotides:
        EIIPxyz[triN] = EIIP_dict[triN[0]] + EIIP_dict[triN[1]] + EIIP_dict[triN[2]]

    encodings = []

    for i in fastas:
        sequence = i.strip()
        code = []
        trincleotide_frequency = TriNcleotideComposition(sequence, base)
        code = code + [EIIPxyz[triN] * trincleotide_frequency[triN] for triN in trincleotides]
        encodings.append(code)

    return np.array(encodings)


chemical_property = {
    'A': [1, 1, 1],
    'C': [0, 1, 0],
    'G': [1, 0, 0],
    'T': [0, 0, 1],
    'U': [0, 0, 1],
    '-': [0, 0, 0],
}





def main1(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = train2DNC(read_nucleotide_sequences(pathTest + TestPosName))
    fea = "DNC.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut + TestPosName[:-6] + fea, header=None, index=None)
    print(fea + "_end.....")

def main2(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = PseEIIP(read_fasta(pathTest + TestPosName))
    fea = "PseEIIP.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut +TestPosName[:-6] + fea, header=None, index=None)
    print(fea + "_end.....")

def main3(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = train2CKSNAP(read_nucleotide_sequences(pathTest + TestPosName))
    fea = "CKSNAP.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut + TestPosName[:-6] + fea, header=None, index=None)
    print(fea + "_end.....")

def main4(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = TNC(read_fasta(pathTest + TestPosName))
    fea = "TNC.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut + TestPosName[:-6]+ fea, header=False, index=None)
    print(fea + "_end.....")



def orderPse(name,TestPosName):
    os.system(
        "python3 /home/lijing/EP3/pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA PC-PseDNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"PC-PseDNC.txt -w 0.1 -k 1")
    os.system(
        "python3 /home/lijing/EP3/pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA PC-PseTNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"PC-PseTNC.txt -w 0.1 -k 1")
    os.system(
        "python3 /home/lijing/EP3/pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA SC-PseTNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"SC-PseTNC.txt -w 0.1 -k 1")
    os.system(
        "python3 /home/lijing/EP3/pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA SC-PseDNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"SC-PseDNC.txt -w 0.1 -k 1")

    os.system(
    "python3 /home/lijing/EP3/pseInOne/acc.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA DACC -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"DACC.txt")



import numpy as np
import pandas as pd
def txt2csv(name,TestPosName):
    features = [
        "PC-PseDNC", "PC-PseTNC", "SC-PseDNC", "SC-PseTNC","DACC"]
    for strFea in features:
        pathTest=name+"temp/"
        TesttxtPos = np.loadtxt(pathTest + TestPosName[:-6] + strFea + ".txt")
        header = ["class"]
        for i in range(len(TesttxtPos[0])):
            header.append(str(i))
        classCol = np.array(list(np.ones(len(TesttxtPos),int)))
        # TesttxtPos = np.insert(TesttxtPos, 0, values=classCol, axis=1)
        TesttxtPos = np.insert(TesttxtPos, 0, values=classCol, axis=1)
        TesttxtDF1 = pd.DataFrame(TesttxtPos)
        TesttxtDF1.to_csv(name + "temp/"+TestPosName[:-6] + strFea + ".csv", index=False, encoding='utf-8', header=header)
        print(strFea + "test_end.......")



import sys, os, re
import numpy as np
import  pandas as pd
def TriNcleotideComposition(sequence, base):
    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    tnc_dict = {}
    for triN in trincleotides:
        tnc_dict[triN] = 0
    for i in range(len(sequence) - 2):
        tnc_dict[sequence[i:i + 3]] += 1
    for key in tnc_dict:
       tnc_dict[key] /= (len(sequence) - 2)
    return tnc_dict





import re
import sys, os
import numpy as np
import pandas as pd
import itertools
from collections import Counter
import numpy as np
def TNC(fastas):
    AA = 'ACGT'
    encodings = []

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i


    for i in fastas:
        sequence = i.strip()
        code = []
        tmpCode = [0] * 64
        for j in range(len(sequence) - 3 + 1):
            tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] = \
                tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] +1

        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return np.array(encodings)



def mainFeature(TestPosName):
    name = "/home/lijing/SubLocEP/"
    ############################ Data
    ####################### featurextract
    main1(name, TestPosName)
    main2(name, TestPosName)
    main3(name, TestPosName)
    main4(name, TestPosName)
    orderPse(name, TestPosName)
    txt2csv(name, TestPosName)

def plk(file):
    df= pd.read_pickle(file)
    return df

def matrix1(file):
    my_matrix = np.loadtxt(open(file,"rb"),delimiter=",",skiprows=0)
    return my_matrix


def matrix2(file):
    my_matrix = np.loadtxt(open(file, "rb"), delimiter=",", skiprows=1)
    X_train = np.delete(my_matrix, 0, axis=1)
    return X_train

def mainModel1(fea,file_name):
    file_name=file_name[:-6]
    path_model="/home/lijing/SubLocEP/model/"
    path_file="/home/lijing/SubLocEP/temp/"
    model=plk(path_model+fea+".plk")
    X_test=matrix2(path_file+file_name+fea+".csv")
    y_pred=model.predict(X_test)
    y_prob = model.predict_proba(X_test)
    np.savetxt("/home/lijing/SubLocEP/temp/"+file_name+"y_pred" + fea + '.csv', y_pred, delimiter=',')
    np.savetxt("/home/lijing/SubLocEP/temp/"+file_name+"y_prob" + fea + '.csv', y_prob, delimiter=',')

def mainModel2(fea,file_name):
    file_name = file_name[:-6]
    path_model="/home/lijing/SubLocEP/model/"
    path_file="/home/lijing/SubLocEP/temp/"
    model=plk(path_model+fea+".plk")
    X_test=matrix1(path_file+file_name+fea+".csv")
    y_pred=model.predict(X_test)
    y_prob = model.predict_proba(X_test)
    np.savetxt("/home/lijing/SubLocEP/temp/"+file_name+"y_pred" + fea + '.csv', y_pred, delimiter=',')
    np.savetxt("/home/lijing/SubLocEP/temp/"+file_name+"y_prob" + fea + '.csv', y_prob, delimiter=',')

import numpy as np
import sklearn.metrics as metrics
from sklearn.metrics import accuracy_score
def matrix(file):

    my_matrix = np.loadtxt(open(file,"rb"),delimiter=",",skiprows=0)
    return my_matrix




def category(av):
    confidence=np.max(av,axis=1)
    index = list(np.argmax(av, axis=1))
    cate=[]
    for inde in index:
        if inde==0:
            cate.append("Cytoplasm")
        if inde==1:
            cate.append("Endoplasmic Reticulum")
        if inde==2:
            cate.append("Extracellular Region")
        if inde==3:
            cate.append("Mitochondria")
        if inde==4:
            cate.append("Nucleus")
    return cate,confidence



def ave(path):
    bc1=matrix(path+"PC-PseDNC.csv")
    bc2= matrix(path + "PC-PseTNC.csv")
    bc3 = matrix(path + "SC-PseDNC.csv")
    bc4 = matrix(path + "SC-PseTNC.csv")
    bc5 = matrix(path + "DACC.csv")

    sb1=matrix(path + "PseEIIP.csv")
    sb2 = matrix(path + "TNC.csv")
    sb3 = matrix(path + "DNC.csv")
    sb4 = matrix(path + "CKSNAP.csv")

    av1=(bc1+bc2+bc3+bc4+bc5)/5
    av2=(sb1+sb2+sb3+sb4)/5
    cate_av1,confidence_av1=category(av1)
    cate_av2,confidence_av2=category(av2)
    av = (av1 + av2 * 2) / 3
    cate_av, confidence_av= category(av)
    return cate_av1,confidence_av1,cate_av2,confidence_av2,cate_av, confidence_av



def single(path):
    bc1 = (matrix(path + "PC-PseDNC.csv"))
    bc2 =  (matrix(path + "PC-PseTNC.csv"))
    bc3 =  (matrix(path + "SC-PseDNC.csv"))
    bc4 =  (matrix(path + "SC-PseTNC.csv"))
    bc5 =  (matrix(path + "DACC.csv"))

    sb1 =  (matrix(path + "PseEIIP.csv"))
    sb2 =  (matrix(path + "TNC.csv"))
    sb3 =  (matrix(path + "DNC.csv"))
    sb4 =  (matrix(path + "CKSNAP.csv"))
    cate_sb1,confidence_sb1=category(sb1)
    cate_sb2, confidence_sb2 = category(sb2)
    cate_sb3, confidence_sb3 = category(sb3)
    cate_sb4, confidence_sb4 = category(sb4)

    cate_bc1, confidence_bc1 = category(bc1)
    cate_bc2, confidence_bc2 = category(bc2)
    cate_bc3, confidence_bc3 = category(bc3)
    cate_bc4, confidence_bc4 = category(bc4)
    cate_bc5, confidence_bc5 = category(bc5)
    return cate_sb1,confidence_sb1,cate_sb2, confidence_sb2,cate_sb3, confidence_sb3,cate_sb4, confidence_sb4,cate_bc1, \
           confidence_bc1,cate_bc2, confidence_bc2,cate_bc3, confidence_bc3,cate_bc4, confidence_bc4,cate_bc5, confidence_bc5


def final(path,filename,seq_name):
    cate_sb1, confidence_sb1, cate_sb2, confidence_sb2, cate_sb3, confidence_sb3, cate_sb4, confidence_sb4, cate_bc1, \
    confidence_bc1, cate_bc2, confidence_bc2, cate_bc3, confidence_bc3, cate_bc4, confidence_bc4, cate_bc5, confidence_bc5=single(path+filename+ "y_prob")
    cate_av1, confidence_av1, cate_av2, confidence_av2,cate_av, confidence_av=ave(path+filename+ "y_prob")
    matrixTemp_sb1 = np.vstack((np.array(cate_sb1), np.array(confidence_sb1)))
    matrixTemp_sb2 = np.vstack((np.array(cate_sb2), np.array(confidence_sb2)))
    matrixTemp_sb3 = np.vstack((np.array(cate_sb3), np.array(confidence_sb3)))
    matrixTemp_sb4 = np.vstack((np.array(cate_sb4), np.array(confidence_sb4)))
    matrixTemp_bc1 = np.vstack((np.array(cate_bc1), np.array(confidence_bc1)))
    matrixTemp_bc2 = np.vstack((np.array(cate_bc2), np.array(confidence_bc2)))
    matrixTemp_bc3 = np.vstack((np.array(cate_bc3), np.array(confidence_bc3)))
    matrixTemp_bc4 = np.vstack((np.array(cate_bc4), np.array(confidence_bc4)))
    matrixTemp_bc5 = np.vstack((np.array(cate_bc5), np.array(confidence_bc5)))
    matrixTemp_av1 = np.vstack((np.array(cate_av1), np.array(confidence_av1)))
    matrixTemp_av2 = np.vstack((np.array(cate_av2), np.array(confidence_av2)))
    matrixTemp_av = np.vstack((np.array(cate_av), np.array(confidence_av)))
    matrixTemp_sb12 = np.vstack((np.array(matrixTemp_sb1), np.array(matrixTemp_sb2)))
    matrixTemp_sb34 = np.vstack((np.array(matrixTemp_sb3), np.array(matrixTemp_sb4)))
    matrixTemp_bc12 = np.vstack((np.array(matrixTemp_bc1), np.array(matrixTemp_bc2)))
    matrixTemp_bc34 = np.vstack((np.array(matrixTemp_bc3), np.array(matrixTemp_bc4)))
    matrixTemp_av12 = np.vstack((np.array(matrixTemp_av1), np.array(matrixTemp_av2)))
    matrixTemp_sb1234 = np.vstack((np.array(matrixTemp_sb12), np.array(matrixTemp_sb34)))
    matrixTemp_bc1234 = np.vstack((np.array(matrixTemp_bc12), np.array(matrixTemp_bc34)))
    matrixTemp_sb1234_bc1234 = np.vstack((np.array(matrixTemp_sb1234), np.array(matrixTemp_bc1234)))
    matrixTemp_sb1234_bc12345 = np.vstack((np.array(matrixTemp_sb1234_bc1234), np.array(matrixTemp_bc5)))
    matrixTemp_sb1234_bc12345_av12 = np.vstack((np.array(matrixTemp_sb1234_bc12345), np.array(matrixTemp_av12)))
    matrixTemp_sb1234_bc12345_av12_av = np.vstack((np.array(matrixTemp_sb1234_bc12345_av12), np.array(matrixTemp_av)))
    header=["seq_name","PseEIIP_LightGBM_label","PseEIIP_LightGBM_score","TNC_LightGBM_label","TNC_LightGBM_score",
            "DNC_LightGBM_label", "DNC_LightGBM_score","CKSNAP_LightGBM_label", "CKSNAP_LightGBM_score",
            "PCPseDNC_LightGBM_label","PCPseDNC_LightGBM_score",
            "PCPseTNC_LightGBM_label", "PCPseTNC_LightGBM_score","SCPseDNC_LightGBM_label","SCPseDNC_LightGBM_score",
            "SCPseTNC_LightGBM_label", "SCPseTNC_LightGBM_score","DACC_LightGBM_label","DACC_LightGBM_score",
            "Sequence-based_single-layer_label", "Sequence-based_single-layer_LightGBM_score",
            "Physicochemical property_single-layer_label","Physicochemical property_single-layer_LightGBM_score",
            "two-layer_model_label", "two-layer_model_score"]
    matrixTemp = np.vstack((seq_name, matrixTemp_sb1234_bc12345_av12_av))
    matrixTemp=matrixTemp .T
    b=pd.DataFrame(matrixTemp)
    b.to_csv("/home/lijing/SubLocEP/static/Result/" +filename + '.fasta_results.csv',index=False,header=header)


def get_name(file):
    fr=open("/home/lijing/SubLocEP/Data/"+file,"r")
    lines=fr.readlines()
    name=[]
    for line in lines:
        if ">" in line:
            name.append(line.replace("\n",""))
    return name


def main(file_name):
    mainFeature(file_name)
    feas = ["DNC", "CKSNAP", "DACC",
            "PC-PseDNC", "PC-PseTNC", "SC-PseDNC", "SC-PseTNC"]
    for fea in feas:
        mainModel1(fea, file_name)
        print(fea)
    feas2 = ["TNC", "PseEIIP", ]
    for fea2 in feas2:
        mainModel2(fea2, file_name)
        print(fea2)
    name = get_name(file_name)
    print(len(name))
    final(path="/home/lijing/SubLocEP/temp/", filename=file_name[:-6] , seq_name=name)




if __name__ == '__main__':
    
    # main("cytoplasm_all_tool_cdhit.txt")
    # main("cytosol_all_tool_cdhit.txt")
    
    # main("ER_all_tool_cdhit.txt")
    # main("exosome_all_tool_cdhit.txt")
    main("Extracellular_region_all_tool_cdhit.txt")
    main("Mitochondria_all_tool_cdhit.txt")
    main("Nucleus_all_tool_cdhit.txt")
    main("posterior_all_tool_cdhit.txt")
    main("Pseudopodium_all_tool_cdhit.txt")
    main("Ribosome_all_tool_cdhit.txt")
