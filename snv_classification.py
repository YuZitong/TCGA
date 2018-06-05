# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-02-12 22:04:18
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-02-12 23:15:23
'''
# pylint: disable=invalid-name

from __future__ import division
import sys
import os
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from sklearn import cross_validation
from snv_pass_site_advanced import workflow

dis_list = ['BLCA', 'BRCA', 'CESC', 'COAD',
            'GBM', 'HNSC', 'KIRC', 'LGG',
            'LIHC', 'LUAD', 'LUSC', 'OV',
            'PRAD', 'SKCM', 'STAD', 'THCA', 'UCEC']

pass_site = []
with open('snv_pass_site_100000_advanced') as ps:
    i = 0
    for line in ps:
        pass_site.append('\t'.join(line.split(',')[0:2]))
        i += 1
        sys.stdout.write('\r' + str(i))
print '\npass site got.'

vcf_file = []
os.system('less snv_pass_site_match_100000_advanced.csv | head -n 10000 > vcf_files')
with open('vcf_files') as vcf:
    i = 0
    for line in vcf:
        vcf_file.append(line.rstrip('\n'))
        i += 1
        sys.stdout.write('\r' + str(i))
print '\nvcf files got.'
os.system('rm vcf_files')

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

print 'start getting datas and labels...'
print 'It may take few time. Please wait...'
d = []
l = []
i = 0
for item in vcf_file:
    wf = workflow(item.split(',')[0])
    with open(item.split(',')[0]) as sample:
        for count, li in enumerate(open(item.split(',')[0], 'rU')):
            pass
        dd = [0 for t in xrange(100000)]
        c = 0
        for line in sample:
            if site(line[0:3] == 'chr', wf, line):
                try:
                    index = pass_site.index('\t'.join(line.split('\t')[0:2]))
                    dd[index] = 1
                    c += 1
                    sys.stdout.write('\r[' + int(float(c/count) * 30) * '#' + int(float((count - c)/count) * 30) * ' ' + ']')
                except ValueError:
                    c += 1
                    sys.stdout.write('\r[' + int(float(c/count) * 30) * '#' + int(float((count - c)/count) * 30) * ' ' + ']')
                    continue
            else:
                c += 1
                sys.stdout.write('\r[' + int(float(c/count) * 30) * '#' + int(float((count - c)/count) * 30) * ' ' + ']')
                continue
        ll = [0 for y in xrange(17)]
        ll[dis_list.index(item.split('/')[4])] = 1
    d.append(dd)
    l.append(ll)
    i += 1
    sys.stdout.write(str(round(i/10000*100, 4)) + '%')
data = np.array(d)
label = np.array(l)
print '\ncomplete!'

data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data, label,
                                                                                     test_size=0.1,
                                                                                     random_state=0)

# model
model = Sequential()
model.add(Dense(10000, activation='relu', input_shape=(100000,)))
# Dropout  to prevent neural networks from overfitting
model.add(Dropout(0.2))
model.add(Dense(5000, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(17, activation='softmax'))

model.summary()
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy',
              optimizer=sgd,
              metrics=['accuracy'])

model.fit(data_train, labels_train,
          epochs=20,
          batch_size=128)

score = model.evaluate(data_test, labels_test, batch_size=128)
print('Test loss:', score[0])
print('Test accuracy:', score[1])