# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-01-15 21:28:47
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-01-18 21:59:32
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
data_src = './data/DNA_meth_report/'
clinical_src = './data/Clinical/'

os.system('find '+ data_src + '|grep .csv >data_list')
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
data_csv = set()
with open('data_list') as rplist:
    for line in rplist:
        data_csv.add(line.rstrip('\n'))
clinical_xml = {} # e.g. clinical_xml[TCGA-ZJ-AAXB] = *clinical xml path*
with open('clinical_list') as rplist:
    for line in rplist:
        clinical_xml[line.split('.')[-2]] = line.rstrip('\n')

print 'extract data...'
d = []
stg = []
count = 0
all_num = len(data_csv)
for paths in data_csv:
    dd = []
    patient_id = paths.split('_')[-2]
    stgstg = [0 for i in xrange(24)]
    try:
        tree = etree.parse(clinical_xml[patient_id])
        root = tree.getroot()
        print root.nsmap['shared_stage'] # for debuging purpose
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
            continue
        stg.append(stgstg)
        with open(paths) as cc:
            for bv in cc:
                dd.append(float(bv.rstrip('\n')))
        d.append(dd)
        sys.stdout.write(str(stgstg) + '\n') # for debuging purpose
    except IOError:
        print '==================' + 'cannot find patient_id:' + patient_id + '==================='
    except KeyError:
        print '==================' + 'cannot find patient_id:' + patient_id + '==================='
    except ValueError:
        print '==================' + 'no stage information' + '==================='
    except IndexError:
        print '==================' + 'no stage information' + '==================='
    count += 1
    sys.stdout.write('\r'+str(round(count/float(all_num)*100, 2))+'%')
meth_bv = numpy.array(d)
stage = numpy.array(stg)
print '\ncomplete!'

# use 10-fold to apply cross validation
data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(meth_bv,
                                                                                     stage,
                                                                                     test_size=0.1,
                                                                                     random_state=0)

print str(len(data_train))+' train samples'
print str(len(data_test))+' test samples'

# model
model = Sequential()
# the first lay: Dense(5000) is a fully-connected layer with 5000 hidden units
model.add(Dense(5000, input_shape=(25978,)))
model.add(advanced_activations.LeakyReLU(alpha=0.01)) # activation: leakyReLU
# Dropout  to prevent neural networks from overfitting
model.add(Dropout(0.3))
# the second lay: Dense(2000) is a fully-connected layer with 2000 hidden units
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
          batch_size=128)

model.save_weights('stage_weights.h5')

score = model.evaluate(data_test, labels_test, batch_size=128)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
