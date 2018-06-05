# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-02-10 17:35:58
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-02-10 19:30:32
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

from __future__ import division
from collections import OrderedDict
import sys

def get_snv_path(snv_src):
    '''
    get all snv path (.vcf)
    input: snv files' path
    output: a list fulled with snv path (.vcf)
    '''
    snv_list = []
    with open(snv_src) as snv:
        for line in snv:
            snv_list.append(line.rstrip('\n'))
    print 'part 1 finished.'
    return snv_list

def get_chrom_pos_dict(file_list):
    '''
    make a dictionary that match chromosome position and the appear numbers of them
    input: a list of snv path (.vcf)
    output: a chrom_pos dictionary. e.g. reault_dict[chrom1,12345] = 1
    '''
    result_dict = {}
    cmplt = 0
    allsite = len(file_list)
    for vcf in file_list:
        wf = workflow(vcf)
        pick(vcf, result_dict, wf)
        cmplt += 1
        print str(round(cmplt/allsite*100, 2)) + '%'
    print '\npart 2 finished.'
    return result_dict

def workflow(vcf_files):
    with open(vcf_files) as v:
        for line in v:
            if 'fileformat' in line:
                continue
            elif 'fileDate' in line:
                #SomaticSniper, ssc and ss
                return 0
            elif 'ID=indelError' in line:
                #VarScan2, ssc and ss'''
                return 1
            elif 'alt_allele_in_normal' in line:
                #MuTect2, PASS'''
                return 2
            elif 'ID=PASS' in line:
                #MuSE, PASS'''
                return 3

def pick(vcf_files, dic, flag):
    if flag == 0:
        with open(vcf_files) as v:
            for line in v:
                try:
                    ss = line.split('\t')[-1].rstrip('\n').split(':')[-2]
                    ssc = int(line.split('\t')[-1].rstrip('\n').split(':')[-1])
                    if ssc > 14 and ss == '2':
                        try:
                            dic[','.join(line.split('\t')[0:2])] += 1
                        except KeyError:
                            dic[','.join(line.split('\t')[0:2])] = 1
                except:
                    continue
    elif flag == 1:
        with open(vcf_files) as v:
            for line in v:
                try:
                    ss = line.split('\t')[7].split(';')[-2]
                    ssc = int(line.split('\t')[7].split('=')[-1])
                    if ssc > 14 and ss == 'SS=2':
                        try:
                            dic[','.join(line.split('\t')[0:2])] += 1
                        except KeyError:
                            dic[','.join(line.split('\t')[0:2])] = 1
                except:
                    continue
    elif flag == 2:
        with open(vcf_files) as v:
            for line in v:
                if '\tPASS' in line:
                    try:
                        dic[','.join(line.split('\t')[0:2])] += 1
                    except KeyError:
                        dic[','.join(line.split('\t')[0:2])] = 1
    elif flag == 3:
        with open(vcf_files) as v:
            for line in v:
                if '\tPASS\tSOMATIC' in line:
                    try:
                        dic[','.join(line.split('\t')[0:2])] += 1
                    except KeyError:
                        dic[','.join(line.split('\t')[0:2])] = 1
    else:
        print 'ERROR: wrong flag!'
        sys.exit()

def write_csv(input_dict, output_path):
    '''
    write report in csv format
    input: input dictionary, output path
    output: none
    '''
    with open(output_path, 'w') as report:
        sorted_chrom_pos_dic = OrderedDict(sorted(input_dict.items(), key=lambda t: t[1], reverse=True))
        for key, value in sorted_chrom_pos_dic.items():
            report.write(key + ',' + str(value) + '\n')

if __name__ == '__main__':
    src = 'snv_list'
    snv_files = get_snv_path(src)
    chrom_pos_dic = get_chrom_pos_dict(snv_files)
    write_csv(chrom_pos_dic, 'snv_pass_site_advanced.csv')
