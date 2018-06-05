# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-03-22 10:49:08
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-03-22 10:49:08
 *
 * tumours classification workflow using DNA methylation data.
 * input: ./data_src/disease_types/uuid/*.txt
 *        downloaded from TCGA and categorized according to disease types
 * output: 1. pretreated data in ./pretreated_src/disease_types/. Format: .csv
 *         2. evaluation of model( accuracy and loss)
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import os
import glob
import sys
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from sklearn import cross_validation
import numpy as np

def get_disease_type_list(src):
    '''
    return disease type list
    input: data source path
    output: disease type list
    '''
    dis_list = glob.glob(src + '/*')
    for i in xrange(len(dis_list)):
        dis_list[i] = dis_list[i].lstrip('/data/DNA_meth/')
    return dis_list

def get_dis_pat_dic(dislist):
    '''
    * return dictionary between disease type(value) and data path(key)
    * make pretreated data dir
    input: disease type list
    output: dis_pat_dic dictionary
    '''
    dis_pat_dic = {}
    for dis in dislist:
        os.system('mkdir '+ pretreated_src +'/'+ dis)
        txt_list = glob.glob(data_src + '/' + dis +'/*/*8.txt')
        for lst in txt_list:
            if int(lst.split('-')[-4][0:2]) < 10:
                dis_pat_dic[lst.rstrip('\n')] = dis
            elif int(lst.split('-')[-4][0:2]) > 9 and int(lst.split('-')[-4][0:2]) < 20:
                dis_pat_dic[lst.rstrip('\n')] = 'Normal'
    return dis_pat_dic

def get_cpg_site(cpgsite):
    '''
    return cpg site set
    input: cpg list file
    output: cpg_set
    '''
    cpgset = set()
    with open(cpgsite) as cpg:
        for line in cpg:
            cpgset.add(line.rstrip('\n'))
    return cpgset

def cpg_input_files(DisPatDic, CpgSet):
    '''
    write pretreated data
    input: (dis_pat_dic, cpgset)
    output: pretreated data files and their path
    '''
    f_no = 0
    num = {}
    total = len(DisPatDic)
    pre_data = []
    for key, value in DisPatDic.items():
        try:
            num[value] += 1
        except KeyError:
            num[value] = 0
        with open(key) as per:
            cpg_meth_bv = []
            for line in per:
                if line.find('Beta_value') != -1:
                    continue
                else:
                    cpgid = line.rstrip('\n').split('\t')[0]
                    bv = line.rstrip('\n').split('\t')[1]
                    if bv == 'NA':
                        bv = '00.0'
                    if cpgid in CpgSet:
                        cpg_meth_bv.append(bv)
            datapath = pretreated_src + value + '/' + value + '_report_'+ '-'.join(key.split('.')[-3].split('-')[0:3]) + '_' + str(num[value]) +'.csv'
            pre_data.append(datapath)
            with open(datapath, 'w') as report:
                bvs = '\n'.join(cpg_meth_bv)
                report.write(bvs + '\n')
            f_no += 1
            sys.stdout.write('\r[' + int(float(f_no/total) * 150) * '=' + int(float((total - f_no)/total) * 150) * ' ' + ']' + str(round(float(f_no/total)*100, 2)) + '%')
    return pre_data

def get_bv(src):
    '''
    read beta value from ./data/DNA_meth_pretreated/*
    input: .csv files(e.g. BRCA_report_TCGA-BH-A0BQ_452.csv)
    output: a list of beta value
    '''
    with open(src) as f:
        bv = f.readlines()
    for i, v in enumerate(bv):
        bv[i] = float(v.rstrip('\n'))
    return bv

def get_label(src):
    '''
    determine tumor type
    input: .csv file's name (e.g. BRCA_report_TCGA-BH-A0BQ_452.csv)
    output: a one-hot-key list represents different tumor type
    '''
    t = src.split('/')[-2]
    label = [0 for i in xrange(18)]
    label[disease_list.index(t)] = 1
    return label

def input_data(file_list):
    '''
    get input data for network training
    input: file list
    output: an array of data for training
    '''
    print 'Preparing input data...'
    print 'It may takes little time. Please be patient.'
    dd = []
    i = 0
    num = len(file_list)
    for item in file_list:
        dd.append(get_bv(item))
        i += 1
        sys.stdout.write('\r[' + int(float(i/num) * 150) * '#' + int(float((num - i)/num) * 150) * ' ' + ']' + str(round(float(i/num)*100, 2)) + '%')
    d = np.array(dd)
    return d

def input_label(file_list):
    '''
    get input label for network training
    input: file list
    output: an array of labels for training
    '''
    print 'Preparing input labels...'
    print 'It may takes little time. Please be patient.'
    ll = []
    i = 0
    num = len(file_list)
    for item in file_list:
        ll.append(get_label(item))
        i += 1
        sys.stdout.write('\r[' + int(float(i/num) * 150) * '#' + int(float((num - i)/num) * 150) * ' ' + ']' + str(round(float(i/num)*100, 2)) + '%')
    l = np.array(ll)
    return l

def build_model(dim0, dim1, dim2, dim3, dropout_rate):
    '''
    build and compile network model
    input:
        dim0: input shape
        dim1: first layer unit
        dim2: second layer unit
        dim3: label number
        dropout rate: drop dropout_rate% input units
    output: model
    '''
    model = Sequential()
    model.add(Dense(dim1, activation='relu', input_shape=(dim0,)))
    # Dropout  to prevent neural networks from overfitting
    model.add(Dropout(dropout_rate))
    model.add(Dense(dim2, activation='relu'))
    model.add(Dropout(dropout_rate))
    model.add(Dense(dim3, activation='softmax'))
    model.summary()
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy',
                  optimizer=sgd,
                  metrics=['accuracy'])
    return model

if __name__ == '__main__':

    data_src = './data/DNA_meth/'
    pretreated_src = './data/DNA_meth_pretreated_test/'
    os.system('mkdir ' + pretreated_src)

    # get disease types according to the directories under data_src
    disease_list = get_disease_type_list(data_src)
    disease_list.append('Normal')

    # matching file with its disease type
    disease_path_dic = get_dis_pat_dic(disease_list)

    # load cpg sites
    cpg_set = get_cpg_site('cpg_list')

    print 'Pretreating input data...'
    data_csv_list = cpg_input_files(disease_path_dic, cpg_set)
    print 'Input date have been pretreated! Path:' + pretreated_src

    # get trainging data and labels
    labels = input_label(data_csv_list)
    data = input_data(data_csv_list)

    # shuffle
    index = np.arange(len(labels))
    np.random.shuffle(index)
    data = data[index]
    labels = labels[index]

    # use 10-fold to apply cross validation
    data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data, labels,
                                                                                         test_size=0.1,
                                                                                         random_state=0)
    print str(len(data_train))+'train samples'
    print str(len(data_test))+'test samples'

    # training
    meth_model = build_model(25978, 2000, 128, 18, 0.2)
    meth_model.fit(data_train, labels_train,
                   epochs=20,
                   batch_size=128)

    # evaluation
    score = meth_model.evaluate(data_test, labels_test, batch_size=128)
    print 'Test loss:', score[0]
    print 'Test accuracy:', score[1]
