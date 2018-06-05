# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-05-30 17:36:30
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-05-30 17:36:30
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
import plotly.plotly as py
import plotly.graph_objs as go
import pandas as pd

dis_list = ['BLCA', 'BRCA', 'CESC', 'COAD', 'GBM', 'HNSC', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PRAD', 'SKCM', 'STAD', 'THCA', 'UCEC', 'Normal']

z = []

true_label = pd.read_csv('label_test.csv', header=None)
predict_label = pd.read_csv('prediction.csv', header=None)

match = {}
for index in true_label.index:
    try:
        match[true_label.loc[index].values.tolist()[0]].append(predict_label.loc[index].values.tolist()[0])
    except KeyError:
        match[true_label.loc[index].values.tolist()[0]] = [predict_label.loc[index].values.tolist()[0]]

for value in match.values():
    correct = [0 for i in xrange(18)]
    for pr in value:
        correct[pr] = correct[pr] + 1.0/len(value)
    z.append(list(correct))

data = [
    go.Heatmap(
        z=z,
        x=dis_list,
        y=dis_list,
        colorscale='Viridis',
    )
]

layout = go.Layout(
    title='Detail of Prediction',
    xaxis = dict(ticks='', nticks=36),
    yaxis = dict(ticks='' )
)

fig = go.Figure(data=data, layout=layout)
py.iplot(fig, filename='prediction-heatmap')