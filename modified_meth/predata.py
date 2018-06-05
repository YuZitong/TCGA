import os
import pandas as pd

data_src = './DNA_meth'
report_src = './data_test'
os.system('mkdir ' + report_src)
dis_string = ['BLCA', 'BRCA', 'CESC', 'COAD', 'GBM', 'HNSC', 'KIRC', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PRAD', 'SKCM', 'STAD', 'THCA', 'UCEC', 'Normal']
blist = pd.read_csv('./cpg_list', header=None)
#os.system('find ' + data_src  + ' -name "*8.txt" > now_list')
clist = []
count = 0
i = 0

with open('now_list') as lis:
    for lst in lis:
        if int(lst.split('-')[-4][0:2]) == 1:
            label = [0] * 18
            label[dis_string.index(lst.split('/')[2])] = 1

            clist = pd.read_table('./' + lst.rstrip('\n'), usecols=['Composite Element REF', 'Beta_value'])
            beta = []
            new_label = []
            new_label.append(label)

            index_deta = clist['Composite Element REF'].isin(blist[0])
            beta = clist[index_deta]

            if count == 0:
                test = pd.DataFrame(data=beta.loc[:, 'Beta_value'])
                test1 = pd.DataFrame(data=new_label)
                count = 1
            else:
                test[count] = beta.loc[:, 'Beta_value']
                test1 = pd.concat([test1, pd.DataFrame(new_label)])
                count = count + 1
        elif int(lst.split('-')[-4][0:2]) >= 10 and int(lst.split('-')[-4][0:2]) <= 19:
            label = [0] * 18
            label[17] = 1

            clist = pd.read_table('./' + lst.rstrip('\n'), usecols=['Composite Element REF', 'Beta_value'])
            beta = []
            new_label = []
            new_label.append(label)

            index_deta = clist['Composite Element REF'].isin(blist[0])
            beta = clist[index_deta]

            if count == 0:
                test = pd.DataFrame(data=beta.loc[:, 'Beta_value'])
                test1 = pd.DataFrame(data=new_label)
                count = 1
            else:
                test[count] = beta.loc[:, 'Beta_value']
                test1 = pd.concat([test1, pd.DataFrame(new_label)])
                count = count + 1
        print lst
        print i/9960.0
    test.to_csv(report_src + '/data.csv', header=None, index=False)
    test1.to_csv(report_src + '/label.csv', header=None, index=False)















