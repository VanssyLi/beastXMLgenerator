import module1_Taxa as mTaxa
import module2_Seq as mSeq
import module3_treeModel as mTreeModel
import module4_popModel as mPop
import module5_clockModel as mClock
import module6_subsModel as mSubs
import module7_likelihood as mLikeli
import module8_operators as mOperator
import module9_priors as mPriors
import module10_logs as mLog
import functions_v2 as f
import math
import xml.etree.ElementTree as ET


fasta = '/Users/vanssyli/Master/SU/CPG/BEASTdata/Mexican_Mammoths/mammoths_mexican_and_10.align_mafft_23463.fasta'
priors_table = '/Users/vanssyli/Master/SU/CPG/Scripts/beast1XMLgenerater/data/Priors_Maxican_Mammoth_10_23463_logNormal.csv'
gff = '/Users/vanssyli/Master/SU/CPG/Scripts/Liftoff/mammoth_maxican_10.liftoff.gff3'
partition_file = '/Users/vanssyli/Master/SU/CPG/Scripts/beast1XMLgenerater/data/partitions.xlsx'
log_name = 'mammoths_mexican_10_23463_logNormal_10million_NoPartition'

taxon_set = [['MX'],
             ['ND']]
ND_list = ['ND']

partition_list = ['tRNA', 'rRNA', 'CDS', 'D_loop']
split_partition = True


beast_setting = f.Xml(units='years',
                      population_model='skygrid',
                      root_height_lower=1000000.0,
                      root_height_mean=5300000.0,
                      root_height_stdev=500000.0,
                      root_height_offset=0.0)

skygrid_cutoff = (int(list(str(beast_setting.root_height_mean))[0]) + 1) * \
                 (math.pow(10, len(str(round(beast_setting.root_height_mean)))))

chainLength = 10000000
log_every = 1000




# Module 1 -- Introduce all the taxa and taxon sets
beast = ET.Element('beast')
tree = ET.ElementTree(beast)
mTaxa.intro_taxa(beast, fasta, priors_table)
mTaxa.intro_taxonSet(beast, taxon_set, fasta)

# Module 2 -- Introduce aligned sequences with/without partitions
mSeq.build_partitionSeq(beast, partition_file, fasta, gff, split_partition)
mSeq.build_patterns(beast, partition_file, split_partition)

# Module 3 -- Introduce the tree model
mTreeModel.build_initialTree(beast, units=beast_setting.units)
mTreeModel.build_treeModel(beast, fasta)

mTreeModel.calc_tip_branchLen(beast)
mTreeModel.calc_taxonSet_branchLen(beast, taxon_set)

# Module 4 -- Define the population model
mPop.skygrid(beast, mean=beast_setting.root_height_mean)

# Module 5 -- Define the clock model
mClock.strict_clock(beast, partition_file, split_partition)

# Module 6 -- Define the substitution model
mSubs.HKY(beast, partition_file, split_partition)

# Module 7 -- Likelihoods
mLikeli.calc_likelihood(beast, partition_file, split_partition)

# Module 8 -- Operators
mOperator.build_operator(beast, partition_file, fasta, split_partition)

# Define MCMC
MCMC = ET.SubElement(beast, 'mcmc', attrib={'id': 'mcmc', 'chainLength': str(chainLength), 'autoOptimize': 'true'})
joint = ET.SubElement(MCMC, 'joint', attrib={'id': 'joint'})
jointPrior = ET.SubElement(joint, 'prior', attrib={'id': 'prior'})
f.make_comment(beast, ' Define MCMC ')

# Module 9 -- Priors
mPriors.define_operators(jointPrior, joint, MCMC, fasta, priors_table, partition_file, split_partition,
                         rootHeight_mean=beast_setting.root_height_mean,
                         rootHeight_stdev=beast_setting.root_height_stdev,
                         rootHeight_offset=beast_setting.root_height_offset)

# Module 10 -- Write Logs
mLog.to_screen(MCMC, partition_file, log_every, split_partition)
mLog.to_file(MCMC, log_every, log_name, taxon_set, partition_file, fasta, split_partition)
mLog.trees_to_file(MCMC, log_every, log_name, partition_file, split_partition)
mLog.report(beast)

# Test and store the XML file
f.pretty_xml(beast)
# ET.dump(beast)

tree.write('{}.xml'.format(log_name), encoding='utf-8', xml_declaration=True)




