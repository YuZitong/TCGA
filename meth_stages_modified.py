# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-03-08 19:59:46
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-03-20 16:49:32
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import sys
import glob
from keras.models import Sequential
from keras.layers import Dense, Dropout, advanced_activations
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping
from sklearn import cross_validation
import numpy as np
from meth_drop import input_data, get_bv
from snv_stage import clinical

stage_list = ['stage 0',
              'stage i', 'stage ia', 'stage ib',
              'stage ii', 'stage iia', 'stage iib', 'stage iic',
              'stage iii', 'stage iiia', 'stage iiib', 'stage iiic',
              'stage iv', 'stage iva', 'stage ivc', 'stage ivb',
              'i/ii nos',
              'stage x']

def get_rid_of_na(file_list):
    '''
    get rif of not report and -- stage files from data_list
    input: data_list
    output: clean data_list
    '''
    p2s = clinical('meth_clinical.tsv')
    patient_id = set()
    for item in file_list:
        p = item.split('_')[-2]
        if p2s[p] == 'not reported' or p2s[p] == '--' or p2s[p] == 'stage x':
            patient_id.add(item)
    print len(patient_id)
    for p_id in patient_id:
        file_list.remove(p_id)
    return file_list

def get_stage(file_name):
    '''
    get stage information of a meth file from meth_clinical.csv
    input: file name path
    output: stage name(in stage_list)
    '''
    global benign
    global malignant
    patien_id = file_name.split('_')[-2]
    stage = patient2stage[patien_id]
    stage_index = [0 for i in xrange(2)]
    if stage_list.index(stage) >= 0 and stage_list.index(stage) <= 3:
        stage_index[0] = 1
        benign += 1
    else:
        stage_index[1] = 1
        malignant += 1
    '''
    if stage_list.index(stage) == 0:
        stage_index[0] = 1
    elif stage_list.index(stage) > 0 and stage_list.index(stage) <= 3:
        stage_index[1] = 1
    elif stage_list.index(stage) > 3 and stage_list.index(stage) <= 7:
        stage_index[2] = 1
    elif stage_list.index(stage) > 7 and stage_list.index(stage) <= 11:
        stage_index[3] = 1
    elif stage_list.index(stage) > 11 and stage_list.index(stage) <= 15:
        stage_index[4] = 1
    elif stage_list.index(stage) == 16:
        stage_index[5] = 1
    else:
        stage_index[6] = 1
    '''
    return stage_index

def input_labels(file_list):
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
        ll.append(get_stage(item))
        i += 1
        sys.stdout.write('\r[' + int(float(i/num) * 150) * '#' + int(float((num - i)/num) * 150) * ' ' + ']' + str(round(float(i/num)*100, 2)) + '%')

    l = np.array(ll)

    return l

if __name__ == '__main__':

    data_src = './data/DNA_meth_pretreated/'
    data_list = glob.glob(data_src + '/*/*.csv')
    data_list = get_rid_of_na(data_list)

    patient2stage = clinical('meth_clinical.tsv')
    print 'patient2stage dictionary got.'

    benign = 0
    malignant = 0
    labels = input_labels(data_list)
    data = input_data(data_list)

    # use 10-fold to apply cross validation
    data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data, labels,
                                                                                         test_size=0.1,
                                                                                         random_state=0)
    print str(len(data_train)) + ' train samples'
    print str(len(data_test)) + ' test samples'

    print str(benign) + ' stage 0 and stage i samples.'
    print str(malignant) + ' stage ii, stage iii and stage iv samples.'

    model = Sequential()
    model.add(Dense(5000, input_shape=(25978,)))
    model.add(advanced_activations.LeakyReLU(alpha=0.01)) # activation: leakyReLU
    model.add(Dropout(0.3))
    model.add(Dense(2000))
    model.add(advanced_activations.LeakyReLU(alpha=0.01))
    model.add(Dropout(0.3))
    model.add(Dense(128))
    model.add(advanced_activations.LeakyReLU(alpha=0.01))
    model.add(Dropout(0.3))
    model.add(Dense(2, activation='softmax'))

    model.summary()

    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy',
                  optimizer=sgd,
                  metrics=['accuracy'])

    model.fit(data_train, labels_train,
              epochs=20,
              batch_size=256)

    score = model.evaluate(data_test, labels_test, batch_size=128)
    print('Test loss:', score[0])
    print('Test accuracy:', score[1])
