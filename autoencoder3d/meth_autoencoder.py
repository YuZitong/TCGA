# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-05-03 16:01:58
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-05-03 16:01:58
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import sys
import numpy as np
np.random.seed(1224)
from keras.datasets import mnist
from keras.models import Model
from keras.layers import Dense, Input
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
from sklearn import cross_validation

dis_list = ['BLCA', 'BRCA', 'CESC', 'COAD', 'GBM', 'HNSC', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PRAD', 'SKCM', 'STAD', 'THCA', 'UCEC', 'Normal']

def get_bv(src):
    '''
    read beta value from ./data/DNA_meth_pretreated/*
    input: .csv files(e.g. BRCA_report_TCGA-BH-A0BQ_452.csv)
    output: a list of beta value
    '''
    with open(src) as f:
        bv = f.readlines()
    for i, v in enumerate(bv):
        bv[i] = float(v.rstrip('\n')) - 0.5
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

def build_model(dim0, dim1, dim2, dim3, encoding_dim):
    '''
    build and train network model
    input:
        dim0: input shape
        dim1: first layer unit
        dim2: second layer unit
        dim3: label number
    output: model
    '''
    # in order to plot in a 2D figure
    #encoding_dim = 2

    # in order to plot in a 3D figure
    #encoding_dim = 3

    # this is our input placeholder
    input_meth = Input(shape=(dim0,))

    # encoder layers
    encoded = Dense(dim1, activation='relu')(input_meth)
    encoded = Dense(dim2, activation='relu')(encoded)
    encoded = Dense(dim3, activation='relu')(encoded)
    encoder_output = Dense(encoding_dim)(encoded)

    # decoder layers
    decoded = Dense(dim3, activation='relu')(encoder_output)
    decoded = Dense(dim2, activation='relu')(decoded)
    decoded = Dense(dim1, activation='relu')(decoded)
    decoded = Dense(dim0, activation='tanh')(decoded)

    # construct the autoencoder model
    autoencoder = Model(input=input_meth, output=decoded)

    # construct the encoder model for plotting
    encoder = Model(input=input_meth, output=encoder_output)

    # compile autoencoder
    autoencoder.compile(optimizer='adam', loss='mse')

    return autoencoder, encoder

def label2num(lable, target):
    num = []
    for l in lable:
        num.append(list(l).index(target))
    return np.array(num)

if __name__ == '__main__':

    data_src = './data/DNA_meth_pretreated/'
    data_list = glob.glob(data_src + '/*/*.csv')

    labels = input_label(data_list)
    data = input_data(data_list)

    # use 10-fold to apply cross validation
    data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data, labels,
                                                                                         test_size=0.1,
                                                                                         random_state=0)
    print str(len(data_train))+' train samples'
    print str(len(data_test))+' test samples'

    autoencoder, encoder = build_model(25978, 5000, 5000, 5000, 2000)

    # training
    autoencoder.fit(data_train, data_train,
                    epochs=5,
                    batch_size=256,
                    shuffle=True)

    encoded_imgs = encoder.predict(data_test)
    np.savetxt('encoder_output.csv', encoded_imgs[:, :], delimiter=',')
    np.savetxt('color.csv', label2num(labels_test, 1), delimiter=',')

    '''
    # plotting
    print 'plotting'
    encoded_imgs = encoder.predict(data_test)
    plt.scatter(encoded_imgs[:, 0], encoded_imgs[:, 1], c=label2num(labels_test, 1))
    plt.colorbar()
    print 'drawing'
    plt.savefig("autoencoder.jpg")
    print 'done'
    '''
