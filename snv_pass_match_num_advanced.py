# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-02-10 19:50:20
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-02-10 19:50:20
'''
# pylint: disable=invalid-name


from __future__ import division
import sys
from collections import OrderedDict

snv_35000_list = []
with open('snv_pass_site_35000_advanced') as snv:
    for line in snv:
        snv_35000_list.append(line.rstrip('\n'))
chrom_pos_set = set()
cmplt = 0
allsite = len(snv_35000_list)
for line in snv_35000_list:
    chrom_pos_set.add(','.join(line.split(',')[0:2]))
    cmplt += 1
    print str(round(cmplt/allsite*100, 2)) + '%'
print '100000 sites got.'

files_list = []
with open('snv_list') as files:
    for smp in files:
        files_list.append(smp.rstrip('\n'))
print 'files path got.'

sample_pass_dic = {}
cmplt = 0
allfiles = len(files_list)
print 'matching...'
for f in files_list:
    with open(f) as snv:
        sample_pass_dic[f] = 0
        for line in snv:
            if 'PASS' in line and ','.join(line.split('\t')[0:2]) in chrom_pos_set:
                sample_pass_dic[f] += 1
    cmplt += 1
    sys.stdout.write('\r' + str(round(cmplt/allfiles*100, 2)) + '% ' + f)
print 'match completed.'

with open('snv_pass_site_match_35000_advanced.csv', 'w') as report:
    sorted_sample_pass_dic = OrderedDict(sorted(sample_pass_dic.items(), key=lambda t: t[1], reverse=True))
    for key, value in sorted_sample_pass_dic.items():
        report.write(key + ',' + str(value) + '\n')
