# 文件说明
## meth_classification.py:
环境：python2.7
python依赖库：sys、os、keras、sklearn
功能：以DNA甲基化为数据，用深度学习分类肿瘤。数据所在的路径：/datapool/home/yuzt/graduation_design/data/DNA_meth_report. 训练好的模型保存在当前目录的weights.h5文件里。

## cpg_importance.m:
环境：Matlab
功能：读取训练好的模型（weights.h5），提取前1000个重要的CpG位点。输出shape为(17,1000)的cpg_importance.csv，CpG位点重要程度随列数递减。

## meth_stages.py:
环境：python2.7
python依赖库：sys、os、keras、sklearn、lxml
功能：以DNA甲基化为数据，用深度学习判断肿瘤所在分期。数据所在的路径：/datapool/home/yuzt/graduation_design/data/DNA_meth_report/Clinical. 训练好的模型保存在当前目录的stage_weights.h5文件里。

# cpg_list
统一两种不同的甲基化文件所用的CpG位点。

# exp_gene_clinical.py
环境：python2.7
python依赖库：sys、os、keras、sklearn、lxml、numpy
功能：以gene expression为数据，判断肿瘤所在分期。

# exp_meth_input.py
环境：python2.7
python依赖库：sys、itertools
功能：为gene expression和甲基化数据用作对分期进行预测的数据做预处理。输出数据地址：/datapool/home/yuzt/graduation_design/data/meth_genexp_inputdata/

# genexp_meth_stage.py
环境：python2.7
python依赖库：sys、lxml、numpy、keras、sklearn
功能：用gene expression和甲基化数据用作对分期进行预测。训练好的模型保存在当前目录的stage_weights.h5文件里。

# snv_classification.py
环境：python 2.7
python依赖库：sys、os、numpy、keras、sklearn、snv_pass_site_advanced
功能：用snv里的Raw Simple Somatic Mutation数据进行肿瘤分类。训练数据所在路径：/datapool/home/yuzt/graduation_design/data/SNV/Raw_Simple_Somatic_Mutation

# snv_pass_site_advanced.py
环境：python 2.7
python依赖库：sys、OrderedDict
功能：整理所有vcf文件里可用的体细胞变异位点。数据所在路径：/datapool/home/yuzt/graduation_design/data/SNV/Raw_Simple_Somatic_Mutation，输出文件：snv_pass_site_advanced.csv

# snv_pass_match_num_advanced.py
环境：python 2.7
python依赖库：sys、OrderedDict
功能：整理出各个vcf文件里拥有的标准体细胞变异位点匹配的位点数量，并按数量降序排列。数据所在路径：/datapool/home/yuzt/graduation_design/data/SNV/Raw_Simple_Somatic_Mutation，标准体细胞变异位点表为snv_pass_site_35000_advanced和snv_pass_site_100000_advanced

# snv_stage.py
环境：python 2.7
python依赖库：sys、os、numpy、keras、sklearn、snv_pass_site_advanced
功能：用snv里的Raw Simple Somatic Mutation数据进行肿瘤分期预测。训练数据所在路径：/datapool/home/yuzt/graduation_design/data/SNV/Raw_Simple_Somatic_Mutation

# meth_stages_modified.py
环境：python2.7
python依赖库：sys、os、keras、sklearn
功能：以DNA甲基化为数据，用深度学习判断肿瘤所在分期。分期分为两类，分别是stage 0&stage i为一类， 其余为一类。数据所在的路径：/datapool/home/yuzt/graduation_design/data/DNA_meth_report.

# meth_drop.py
环境：python2.7
python依赖库：sys、os、keras、sklearn
功能：以不完整的片段化的DNA甲基化为数据，用深度学习重新训练模型，分类肿瘤。数据所在的路径：/datapool/home/yuzt/graduation_design/data/DNA_meth_report.

# meth_drop_test.py
环境：python2.7
python依赖库：sys、os、keras、sklearn
功能：以不完整的片段化的DNA甲基化为数据，加载原先训练完毕的模型，进行肿瘤分类的预测。数据所在的路径：/datapool/home/yuzt/graduation_design/data/DNA_meth_report.

# meth_classification_workflow
环境：python2.7
python依赖库：sys、os、keras、sklearn
DNA甲基化预测分类流程。其中training_workflow结尾的是使用整个TCGA数据进行训练的流程，predict结尾的是对输入的甲基化数据进行分类预测， improve结尾的是使用新加入的甲基化数据进行模型的进一步训练。

# modified_meth
环境：python2.7
python依赖库：同meth_classification_workflow
剔除除01标签以外的肿瘤样本，并进行训练、预测、生成混淆矩阵的流程。

# 附录
## Keras 文档
https://keras-cn.readthedocs.io/en/latest/
