# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-03-04 16:28:25
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-03-04 16:28:25
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import random
import sys
import glob
from meth_drop import *
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.optimizers import SGD
from keras.callbacks import EarlyStopping
from sklearn import cross_validation
import numpy as np

def ran(num, rate):
    '''
    generate a num length list whose values are 0 and 1. The proportion is depent on rate.
    input:
        num: length of list
        rate: the proportion of 1 in list
    output: result list
    '''
    x = [1 for t in xrange(num)]
    for i, v in enumerate(x):
        if random.random() > rate:
            x[i] = 0
    return x

if __name__ == '__main__':
    data_src = './data/DNA_meth_report/'
    data_list = glob.glob(data_src + '/*/*.csv')

    labels = input_label(data_list)
    data = input_data(data_list)

    num = len(data)
    rate = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    for rr in rate:
        print 'modifing input data...'
        data_m = data.copy()
        for i, v in enumerate(data):
            data_m[i] = data[i] * np.array(ran(25978, rr))
            sys.stdout.write('\r' + str(round(float(i/num)*100, 2)) + '%')

        model = build_model(25978, 2000, 128, 18, 1)
        model.load_weights('weights.h5')

        tt = 0
        for c in data_m[0]:
            if c == 0:
                tt += 1

        score = model.evaluate(data_m, labels, batch_size=128)
        print '0 proportion', str(tt/25978)
        print 'rate:', str(rr)
        print 'Test loss:', score[0]
        print 'Test accuracy:', score[1]
