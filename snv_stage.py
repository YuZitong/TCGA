# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-02-20 13:48:34
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-02-20 13:48:34
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import sys
import os
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping
from sklearn import cross_validation
from snv_pass_site_advanced import workflow

stages = ['stage 0',
          'stage i', 'stage ia', 'stage ib',
          'stage ii', 'stage iia', 'stage iib', 'stage iic',
          'i/ii nos',
          'stage iii', 'stage iiia', 'stage iiib', 'stage iiic',
          'stage iv', 'stage iva', 'stage ivb', 'stage ivc',
          'stage x']

def site(start, flag, site_line):
    '''
    check the input site
    '''
    if start:
        if flag == 0:
            ss = site_line.split('\t')[-1].rstrip('\n').split(':')[-2]
            ssc = int(site_line.split('\t')[-1].rstrip('\n').split(':')[-1])
            if ssc > 14 and ss == '2':
                return True
            else:
                return False
        elif flag == 1:
            ss = site_line.split('\t')[7].split(';')[-2]
            ssc = int(site_line.split('\t')[7].split('=')[-1])
            if ssc > 14 and ss == '2':
                return True
            else:
                return False
        elif flag == 2:
            if '\tPASS' in line:
                return True
            else:
                return False
        elif flag == 3:
            if '\tPASS\tSOMATIC' in line:
                return True
            else:
                return False
    else:
        return False

def sample_sheet(sheet_path):
    '''
    associate snv files with patient id
    input: sample sheet path
    output: a dictionary whose snv files' path and whose values are keys are patient id
    '''
    snv2patient = {}
    with open(sheet_path) as f:
        for line in f:
            if 'File ID' in line:
                continue
            else:
                snv2patient[line.split('\t')[0]] = line.split('\t')[5].split(',')[0]
    return snv2patient

def clinical(clinical_tsv_path):
    '''
    associate patient id with tumor_stage
    input: snv_clinical.tsv path
    output: a dictionary whose keys are patient id and whose values are tumor stage
    '''
    patient2stage = {}
    with open(clinical_tsv_path) as f:
        for line in f:
            if 'case_id' in line:
                continue
            else:
                patient2stage[line.split('\t')[1]] = line.split('\t')[11]
    return patient2stage

def passsite(path):
    '''
    get pass site list
    input: the pass site file(e.g. snv_pass_site_100000_advanced, etc)
    output: a list of pass site
    '''
    pass_site = []
    with open(path) as ps:
        for line in ps:
            pass_site.append('\t'.join(line.split(',')[0:2]))
    print 'pass site got.'
    return pass_site

def data(snv_path, input_dim):
    '''
    get training data of one sample
    input: modified-top-10000 snv_path
    output: network training data array of one sample
    '''
    wf = workflow(snv_path)
    dd = [0 for t in xrange(input_dim)]
    c = 0
    with open(snv_path) as f:
        for count, li in enumerate(open(snv_path, 'rU')):
            pass
        for line in f:
            if site(line[0:3] == 'chr', wf, line):
                try:
                    index = pass_site.index('\t'.join(line.split('\t')[0:2]))
                    dd[index] = 1
                    sys.stdout.write('\r[' + int(float(c/count) * 30) * '#' + int(float((count - c)/count) * 30) * ' ' + ']')
                    c += 1
                except ValueError:
                    sys.stdout.write('\r[' + int(float(c/count) * 30) * '#' + int(float((count - c)/count) * 30) * ' ' + ']')
                    c += 1
                    continue
            else:
                sys.stdout.write('\r[' + int(float(c/count) * 30) * '#' + int(float((count - c)/count) * 30) * ' ' + ']')
                c += 1
                continue
    return dd

def label(snv_path):
    '''
    get label of one sample
    input: modified-top-10000 snv_path
    output: network label array of one sample. if stage is not reported or --, return 0
    '''
    ll = [0 for i in xrange(len(stages))]
    patient_id = snv2patient[snv_path.split('/')[-1].split('.')[-2]]
    stage = patient2stage[patient_id]
    if stage == 'not reported' or stage == '--':
        return 0
    index = stages.index(stage)
    ll[index] = 1
    return ll

if __name__ == '__main__':
    snv2patient = sample_sheet('snv_sample_sheet.tsv')
    print 'snv2patient dictionary got.'
    patient2stage = clinical('snv_clinical.tsv')
    print 'patient2stage dictionary got.'
    pass_site = passsite('snv_pass_site_100000_advanced')
    print str(len(pass_site)) + ' pass sites got.'

    vcf_file = []
    os.system('less snv_pass_site_match_100000_advanced.csv | head -n 10000 > vcf_files')
    with open('vcf_files') as vcf:
        for line in vcf:
            vcf_file.append(line.rstrip('\n'))
    print 'vcf files got.'
    os.system('rm vcf_files')

    print 'start getting datas and labels...'
    print 'It may take few time. Please wait...'
    d = []
    l = []
    j = 0
    for item in vcf_file:
        path = item.split(',')[0]
        flag = label(path)
        if flag == 0:
            j += 1
            continue
        l.append(flag)
        d.append(data(path, 100000))
        j += 1
        sys.stdout.write(str(round(j/10000*100, 4)) + '%')
    fnldata = np.array(d)
    fnllabel = np.array(l)

    data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(fnldata, fnllabel,
                                                                                         test_size=0.1,
                                                                                         random_state=0)

    # model
    model = Sequential()
    model.add(Dense(10000, activation='relu', input_shape=(100000,)))
    # Dropout  to prevent neural networks from overfitting
    model.add(Dropout(0.2))
    model.add(Dense(2000, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(len(stages), activation='softmax'))

    model.summary()
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy',
                  optimizer=sgd,
                  metrics=['accuracy'])

    early_stopping = EarlyStopping(monitor='acc', patience=2)
    model.fit(data_train, labels_train,
              epochs=20,
              batch_size=128,
              callbacks=[early_stopping])

    score = model.evaluate(data_test, labels_test, batch_size=128)
    print('Test loss:', score[0])
    print('Test accuracy:', score[1])
