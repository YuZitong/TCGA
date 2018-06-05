# -*- coding: utf-8 -*-
'''
* @Author: Yu Zitong
* @Date: 2018-01-06 14:22:11
* @Last Modified by:   Yu Zitong
* @Last Modified time: 2018-01-14 12:25:11
'''
# pylint: disable=invalid-name

import sys
import os
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from sklearn import cross_validation

# get meth beta_value data
import numpy as np

data_src = './data/DNA_meth_report/'

os.system('find '+ data_src + '|grep .csv >data_list')

dis = []
with open('./dis_list') as dslst:
    for line in dslst:
        dis.append(line.rstrip('\n'))
dis.append('Normal')

print 'extract data path...'
data_csv = set()
with open('data_list') as rplist:
    for line in rplist:
        data_csv.add(line.rstrip('\n'))

print 'extract data...'
d = []
l = []
count = 0
all_num = len(data_csv)
for paths in data_csv:
    dd = []
    ll = [0 for i in xrange(18)]
    with open(paths) as cc:
        for bv in cc:
            dd.append(float(bv.rstrip('\n')))
    d.append(dd)
    ll[dis.index(paths.split('/')[-2])] = 1
    l.append(ll)
    count += 1
    sys.stdout.write('\r'+str(round(count/float(all_num)*100, 2))+'%')
data = np.array(d)
labels = np.array(l)
print '\ncomplete!'

# use 10-fold to apply cross validation
data_train, data_test, labels_train, labels_test = cross_validation.train_test_split(data, labels,
                                                                                     test_size=0.1,
                                                                                     random_state=0)

print str(len(data_train))+'train samples'
print str(len(data_test))+'test samples'

# model
model = Sequential()
# the first lay: Dense(2000) is a fully-connected layer with 2000 hidden units
model.add(Dense(2000, activation='relu', input_shape=(25978,)))
# Dropout  to prevent neural networks from overfitting
model.add(Dropout(0.2))
# the second lay: Dense(128) is a fully-connected layer with 128 hidden units
model.add(Dense(128, activation='relu'))
# Dropout  to prevent neural networks from overfitting
model.add(Dropout(0.2))
# the last lay: Dense(18) is a fully-connected layer with 17 output units
model.add(Dense(18, activation='softmax'))

model.summary()
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy',
              optimizer=sgd,
              metrics=['accuracy'])

model.fit(data_train, labels_train,
          epochs=20,
          batch_size=128)

model.save_weights('weights.h5')

score = model.evaluate(data_test, labels_test, batch_size=128)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
