# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-03-01 13:49:00
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-03-03 14:31:23
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import sys
import glob
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping
from sklearn import cross_validation
import numpy as np

dis_list = ['BLCA', 'BRCA', 'CESC', 'COAD', 'GBM', 'HNSC', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PRAD', 'SKCM', 'STAD', 'THCA', 'UCEC', 'Normal']

def get_bv(src):
    '''
    read beta value from ./data/DNA_meth_report/*
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

    label[dis_list.index(t)] = 1

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
    model.add(Dense(dim3, activation='softmax'))

    model.summary()
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy',
                  optimizer=sgd,
                  metrics=['accuracy'])

    return model

if __name__ == '__main__':
    data_src = './data/DNA_meth_report/'
    data_list = glob.glob(data_src + '/*/*.csv')

    labels = input_label(data_list)
    data = input_data(data_list)

    # use 10-fold to apply cross validation
    data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data, labels,
                                                                                         test_size=0.1,
                                                                                         random_state=0)
    print str(len(data_train))+'train samples'
    print str(len(data_test))+'test samples'

    dropout_list = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
    finial_accuracy = []

    for rate in dropout_list:
        model = build_model(25978, 2000, 128, 18, rate)

        early_stopping = EarlyStopping(monitor='acc', patience=1)
        model.fit(data_train, labels_train,
                  epochs=20,
                  batch_size=128)

        score = model.evaluate(data_test, labels_test, batch_size=128)
        print('Test loss:', score[0])
        print('Test accuracy:', score[1])
        finial_accuracy.append(score[1])

    print finial_accuracy
