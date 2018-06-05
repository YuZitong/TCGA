# -*- coding: utf-8 -*-
'''
 * @Author: Yu Zitong
 * @Date: 2018-01-23 19:09:01
 * @Last Modified by:   Yu Zitong
 * @Last Modified time: 2018-01-24 11:23:56
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

import sys
import itertools

def file_case():
    '''
    match file names with case ID
    '''
    f2c = {} # f2c[90629e6e-2e1e-4f31-884a-6bca2c24ce11.FPKM.txt.gz] = TCGA-EL-A3ZO
    with open(clinical2sampleID_src) as tsv:
        for line in tsv:
            if 'File ID' in line:
                continue
            elif int(line.split('\t')[-2].split('-')[-1][0:2]) < 10:
                f2c[line.split('\t')[1]] = line.split('\t')[-3]
    return f2c


def get_meth_path(path):
    '''
    get meth file's path
    input: data_list
           e.g. ./data/DNA_meth_report/Normal/Normal_report_TCGA-EL-A3ZO_973.csv
    output: meth_path_dic
            e.g. meth_path_dic[TCGA-EL-A3ZO] = ./data/DNA_meth_report/Normal/Normal_report_TCGA-EL-A3ZO_973.csv
    '''
    output_dic = {}
    with open(path) as meth:
        for line in meth:
            output_dic[line.split('_')[-2]] = line.rstrip('\n')
    return output_dic

def get_gene_exp_path(path):
    '''
    get gene expression file's path
    input: gene_exp_data_list
           e.g. ./data/Exp/Gene_Exp/GBM/eae23ce3-2560-4825-9027-9150b19ad0fd/90629e6e-2e1e-4f31-884a-6bca2c24ce11.FPKM.txt
    output: gene_exp_path_dic
            e.g. gene_exp_path_dic[TCGA-EL-A3ZO] = ./data/Exp/Gene_Exp/GBM/eae23ce3-2560-4825-9027-9150b19ad0fd/90629e6e-2e1e-4f31-884a-6bca2c24ce11.FPKM.txt
    '''
    output_dic = {}
    with open(path) as genexp:
        for line in genexp:
            try:
                output_dic[file2case_dic[line.rstrip('\n').split('/')[-1] + '.gz']] = line.rstrip('\n')
            except KeyError:
                pass
    return output_dic

def match_write(dic1, dic2):
    '''
    match meth with gene exp and write csv
    input: meth_path_dic
           gene_exp_path_dic
    output: (list)content
    '''
    meth_p = ''
    genexp_p = ''
    num = len(dic2)
    cmp = 0
    for key, value in dic2.items():
        if key in dic1.keys():
            meth_p = dic1[key]
            genexp_p = value
            p1 = []
            p2 = []
            with open(meth_p) as part1:
                for line in part1:
                    p1.append(line.rstrip('\n') + ',')
            with open(genexp_p) as part2:
                for line in part2:
                    p2.append(line.rstrip('\n').split('\t')[-1])
            content = map(lambda (a, b): a + b + '\n', itertools.izip_longest(p1, p2, fillvalue = '00.0,'))
            with open('./data/meth_genexp_inputdata/' + key + '.csv', 'w') as output:
                output.writelines(content)
            cmp += 1
            sys.stdout.write('\r' + str(round(float(cmp)/num*100, 2)) + '%')

if __name__ == '__main__':

    clinical2sampleID_src = './data/Exp/Gene_Exp/file2sample.tsv'

    file2case_dic = file_case()

    print 'getting meth file\'s path...'
    meth_path_dic = get_meth_path('./meth_data_list')

    print 'gettings gene expression file\'s path...'
    gene_exp_path_dic = get_gene_exp_path('./gene_exp_data_list')

    print 'matching and writing...'
    match_write(meth_path_dic, gene_exp_path_dic)
