# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-05-28 14:44:51
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-05-28 14:44:51
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import glob
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping
from sklearn import cross_validation
import numpy as np
import pandas as pd

dis_list = ['BLCA', 'BRCA', 'CESC', 'COAD', 'GBM', 'HNSC', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PRAD', 'SKCM', 'STAD', 'THCA', 'UCEC', 'Normal']

def build_model(dim0, dim1, dim2, dim3, dropout_rate):
    '''
    build and compile network model
    input:
        dim0: input shape
        dim1: first layer unit
        dim2: second layer unit
        dim3: label number
        dropout rate: dropout_rate
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
    data_src = './data_test/data.csv'
    label_src = './data_test/label.csv'

    print 'Loading data...'
    data = pd.read_csv(data_src, header=None)
    data = data.T
    data = data.fillna(value=0.0)
    print 'Loading label...'
    label = pd.read_csv(label_src, header=None)

    print 'train-test spliting...'
    # use 10-fold to apply cross validation
    data_train, data_test, label_train, label_test = cross_validation.train_test_split(data, label, test_size=0.1, random_state=0)

    print str(len(data_train))+' train samples'
    print str(len(data_test))+' test samples'

    #model building
    meth_model = build_model(25978, 2000, 128, 18, 0.2)

    early_stopping = EarlyStopping(monitor='acc', patience=3)
    history = meth_model.fit(np.array(data_train), np.array(label_train), epochs=20, batch_size=128, callbacks=[early_stopping])

    meth_predict = meth_model.predict_classes(data_test, batch_size=128, verbose=0)

    labels = []
    for ll in label_test.values:
        labels.append(np.where(ll == 1)[0].tolist()[0])
    labels = np.array(labels)

    np.savetxt('label_test.csv', labels, fmt='%d', delimiter=',', newline='\n')
    np.savetxt('prediction.csv', meth_predict, fmt='%d', delimiter=',', newline='\n')
