# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-01-24 11:24:09
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-01-25 10:29:23
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

import sys
from lxml import etree # for processing XML
import numpy
from keras.models import Sequential
from keras.layers import Dense, Dropout, advanced_activations, Flatten
from keras.optimizers import SGD
from sklearn import cross_validation
from sklearn import preprocessing

def get_data(path):
    '''
    get data
    input: data path
           e.g. ./data/meth_genexp_inputdata/TCGA-GN-A9SD.csv
    output: meth + gene expression data
            e.g. [[00.0,2.34],
                  [0.342,2.34],
                  ... ...
                  [0.435,23.234]]
    '''
    data = numpy.loadtxt(path, dtype=float, delimiter=',')
    return data.tolist()

def get_label(caseid):
    '''
    get labels
    input: case id
           e.g. TCGA-GN-A9SD
    output: stage information
            e.g. [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    '''
    stgstg = [0 for i in xrange(24)]
    try:
        tree = etree.parse(clinical_xml[caseid])
        root = tree.getroot()
        tt = root.nsmap['shared_stage']
        disease = root.findall('admin:admin/admin:disease_code', root.nsmap)
        clnstg = root.findall(str.lower(disease[0].text) +
                              ':patient/shared_stage:stage_event/shared_stage:clinical_stage',
                              root.nsmap)
        pathstg = root.findall(str.lower(disease[0].text) +
                               ':patient/shared_stage:stage_event/shared_stage:pathologic_stage',
                               root.nsmap)
        if clnstg[0].text != None:
            stgstg[find_stage(clnstg[0].text)] = 1
        elif pathstg[0].text != None:
            stgstg[find_stage(pathstg[0].text)] = 1
        else:
            return 'warnings'
    except:
        return 'warnings'
    return stgstg

def find_stage(t):
    '''
    for descending the label's dimension
    input: stage name
    output: index in 'stages'
    '''
    p = 100
    for sub in stages:
        if t in sub:
            p = stages.index(sub)
            break
    return p

def build_train(dtr, ltr, dte, lte):
    '''
    build, train model and evaluate it
    input:
          dtr: train data
          ltr: train labels
          dte: test data
          lte: test labels
    output: model evaluation report and weights.h5
    '''
    model = Sequential()
    model.add(Flatten(input_shape=(60483,2)))
    # the first lay: Dense(5000) is a fully-connected layer with 5000 hidden units
    model.add(Dense(16000))
    model.add(advanced_activations.LeakyReLU(alpha=0.01)) # activation: leakyReLU
    # Dropout  to prevent neural networks from overfitting
    model.add(Dropout(0.3))
    # the second lay: Dense(1000) is a fully-connected layer with 2000 hidden units
    model.add(Dense(2000))
    model.add(advanced_activations.LeakyReLU(alpha=0.01))
    model.add(Dropout(0.3))
    # the third lay: Dense(128) is a fully-connected layer with 128 hidden units
    model.add(Dense(128))
    model.add(advanced_activations.LeakyReLU(alpha=0.01))
    model.add(Dropout(0.3))
    # the last lay: Dense(24) is a fully-connected layer with 24 output units
    model.add(Dense(24, activation='softmax'))

    model.summary()

    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy',
                  optimizer=sgd,
                  metrics=['accuracy'])

    model.fit(dtr, ltr,
              epochs=20,
              batch_size=256)

    model.save_weights('stage_weights.h5')

    score = model.evaluate(dte, lte, batch_size=256)
    print('Test loss:', score[0])
    print('Test accuracy:', score[1])

if __name__ == '__main__':

    clinical_xml = {} # e.g. clinical_xml[TCGA-ZJ-AAXB] = *clinical xml path*
    with open('clinical_list') as rplist:
        for line in rplist:
            clinical_xml[line.split('.')[-2]] = line.rstrip('\n')

    stages = ['Stage I',
              'Stage IA', 'Stage IA1', 'Stage IA2',
              'Stage IB', 'Stage IB1', 'Stage IB2',
              'Stage IC',
              'Stage II',
              'Stage IIA', 'Stage IIA1', 'Stage IIA2',
              'Stage IIB',
              'Stage IIC',
              'Stage III',
              'Stage IIIA',
              'Stage IIIB',
              'Stage IIIC', 'Stage IIIC1', 'Stage IIIC2',
              'Stage IV',
              'Stage IVA',
              'Stage IVB',
              'Stage IVC']

    meth_fpkm = []
    stage = []
    with open('./meth_fpkm_list') as inputlist:
        cmplt = 0
        for line in inputlist:
            if get_label(line.rstrip('.csv\n').split('/')[-1]) != 'warnings':
                meth_fpkm.append(preprocessing.scale(get_data(line.rstrip('\n'))))
                stage.append(get_label(line.rstrip('.csv\n').split('/')[-1]))
            else:
                pass
            cmplt += 1
            sys.stdout.write('No. ' + str(cmplt) + ' :' +str(round(cmplt/80.67, 2)) + '%\n')

    data = numpy.array(meth_fpkm)
    label = numpy.array(stage)

    # use 10-fold to apply cross validation
    data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data,
                                                                                         label,
                                                                                         test_size=0.1,
                                                                                         random_state=0)
    print str(len(data_train))+' train samples'
    print str(len(data_test))+' test samples'

    build_train(data_train, labels_train, data_test, labels_test)
