# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-01-22 20:56:46
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-01-22 20:56:46
'''
# pylint: disable=invalid-name

import os
import sys
from lxml import etree # for processing XML
import numpy
from keras.models import Sequential
from keras.layers import Dense, Dropout, advanced_activations
from keras.optimizers import SGD
from sklearn import cross_validation

# load data
data_src = './data/Exp/Gene_Exp'
clinical_src = './data/Clinical/'

os.system('find '+ data_src + ' -name "*.FPKM.txt" >data_list')
os.system('find ' + clinical_src + ' -name "*.xml" >clinical_list')

## disease list
dis = []
with open('./dis_list') as dslst:
    for line in dslst:
        dis.append(line.rstrip('\n'))

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

print 'extract data path...'
data_txt = set()
with open('data_list') as rplist:
    for line in rplist:
        data_txt.add(line.rstrip('\n'))
clinical_xml = {} # e.g. clinical_xml[TCGA-ZJ-AAXB] = *clinical xml path*
with open('clinical_list') as rplist:
    for line in rplist:
        clinical_xml[line.split('.')[-2]] = line.rstrip('\n')

# match file names with sample ID
file2sample_dic = {}
with open('./data/Exp/Gene_Exp/file2sample.tsv') as tsv:
    for line in tsv:
        if 'File ID' in line:
            continue
        else:
            file2sample_dic[line.split('\t')[1]] = line.split('\t')[-2]

print 'extract data...'
d = []
stg = []
count = 0
all_num = len(data_txt)
label = []
ns = 0
cf = 0
for paths in data_txt:
    if int(file2sample_dic[paths.rstrip('\n').split('/')[-1] + '.gz'].split('-')[-1][0:2]) < 10:
        dd = []
        patient_id = '-'.join(file2sample_dic[paths.rstrip('\n').split('/')[-1] + '.gz'].split('-')[0:3])
        stgstg = [0 for i in xrange(24)]
        try:
            tree = etree.parse(clinical_xml[patient_id])
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
                ns += 1
                count += 1
                continue
            stg.append(stgstg)
            with open(paths) as genexp:
                for fpkm in genexp:
                    dd.append(float(fpkm.rstrip('\n').split('\t')[-1]))
            d.append(dd)
        except KeyError:
            cf += 1
        except IndexError:
            ns += 1
        count += 1
        sys.stdout.write('\r'+str(round(count/float(all_num)*100, 2))+'% ' + str(stgstg))
    else:
        count += 1
exp_gene_fpkm = numpy.array(d)
stage = numpy.array(stg)
print '\ncomplete!' + str(ns) + ' samples: no stage information. ' + str(cf) + ' samples: no clinical records.'

# use 10-fold to apply cross validation
data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(exp_gene_fpkm,
                                                                                     stage,
                                                                                     test_size=0.1,
                                                                                     random_state=0)

print str(len(data_train))+' train samples'
print str(len(data_test))+' test samples'

# model
model = Sequential()
# the first lay: Dense(5000) is a fully-connected layer with 5000 hidden units
model.add(Dense(8000, input_shape=(60483,)))
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

model.fit(data_train, labels_train,
          epochs=20,
          batch_size=256)

model.save_weights('stage_weights.h5')

score = model.evaluate(data_test, labels_test, batch_size=256)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
