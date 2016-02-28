###############################################################################
# Reference implementation based on variable elimination
# author: Jinjing Zhou, Isaac Caswell
# date: Jan 20, 2016

###############################################################################

import copy
import numpy as np
import sys
from factors import *
from binary_tree import BinaryTree
import math
from sys import stdout


def read_nh_file_into_single_line(filename):
    with open(filename) as f:
        str_to_parse = f.read().replace('\n','').replace(' ','').strip(' \t\n\r;')
    return str_to_parse

def read_sequences(filename, animal_to_idx):
    observations = {}
    with open(filename) as f:
        line = f.read()
        raw_data = line.split('>')
        for i in xrange(1, len(raw_data)):
            data = raw_data[i].strip().split('\n')
            observations[animal_to_idx[data[0]]] = ''.join(data[1:])
    return observations

def generate_factor(treeFile, columnFile, prior, cpd, domain):

    # Build the binary tree and find the interior nodes
	tree = BinaryTree()
	tree.parse_from_string(read_nh_file_into_single_line(treeFile))
	interior_nodes = []
	tree.assign_idx(interior_nodes, [0])


	animal_to_idx={}
	idx_to_animal={}
	factors=[]


	def build_prior_factor(idx):
		factor=Factor()
		factor.scope = [idx]
		factor.card = [len(domain)]
		factor.val = prior
		return factor
	def build_evolution_factor(i2,i1):
		factor=Factor()
		factor.scope = [i2,i1]
		factor.card = [len(domain),len(domain)]
		factor.val = cpd
		return factor


	def assign_index(TreeNode):
		animal_to_idx[TreeNode.data]=TreeNode.idx+1
		idx_to_animal[TreeNode.idx+1]=TreeNode.data
		if TreeNode.left!=None:
			f_e=build_evolution_factor(TreeNode.left.idx+1,TreeNode.idx+1)
			f_e.name = TreeNode.left.data+'<--'+TreeNode.data
			factors.append(f_e)
			assign_index(TreeNode.left)
		if TreeNode.right!=None:
			f_e=build_evolution_factor(TreeNode.right.idx+1,TreeNode.idx+1)
			f_e.name = TreeNode.right.data+'<--'+TreeNode.data
			factors.append(f_e)
			assign_index(TreeNode.right)

	f_p=build_prior_factor(tree.idx+1)
	f_p.name=tree.data
	factors.append(f_p)
	assign_index(tree)

	sequences = read_sequences(columnFile,animal_to_idx)
	domain_to_idx = {domain[i]: i for i in xrange(len(domain))}

	return sequences,tree,factors


def tree_order(tree):
	od=[]
	def order(tree):
		if tree.left!=None:
			if tree.left.is_leaf():
				od.append(tree.left.idx+1)
				tree.left=None
			else:
				order(tree.left)
		if tree.right!=None:
			if tree.right.is_leaf():
				od.append(tree.right.idx+1)
				tree.right=None
			else:
				order(tree.right)			
	tree_copy=copy.deepcopy(tree)
	while not tree_copy.is_leaf():
		order(tree_copy)
	od.append(tree_copy.idx+1)
	return od

def VE_log(observation,tree,factors,od):
	for ob in observation.items():
		ch_idx=ob[0]
		for f in factors:
			if ch_idx in f.scope:
				factors.remove(f)
				f=f.observe(ch_idx,domain.index(ob[1]))
				factors.append(f)
				break
	#table=[]
	for ch_idx in od:
		#message=[ch_idx]
		message_factors=[]
		removal=[]
		for f in factors:
			if ch_idx in f.scope:
				if len(f.scope)>1:
					pa_idx=[p for p in f.scope if p!=ch_idx][0]
					#message.append(pa_idx)

				message_factors.append(f)
				removal.append(f)
		for r in removal:
			factors.remove(r)
		#log function here
		
		m=Factor(message_factors[0])
		for i in range(1,len(message_factors)):
			m=m.sum(message_factors[i])

		marginal=copy.deepcopy(m.scope)
		marginal.remove(ch_idx)
		for a in indices_to_assignment(range(np.prod(m.card)),m.card):
			if m.val[tuple(a)]!=0:
				m.val[tuple(a)]=math.exp(m.val[tuple(a)])
		m=m.marginalize_all_but(marginal)
		m.val=np.log(m.val)
		factors.append(m)
		#message.append(m.val)
		#table.append(message)
	#return table,table[-1][-1]
	return m.val


if __name__ == '__main__':
	prior = np.log([0.295, 0.205, 0.205, 0.295])
	cpd = np.log( \
	[[0.831, 0.046, 0.122, 0.053], \
	 [0.032, 0.816, 0.028, 0.076], \
	 [0.085, 0.028, 0.808, 0.029], \
	 [0.052, 0.110, 0.042, 0.842]])
	domain = ['A', 'C', 'G', 'T']   
	sequences,tree,factors=generate_factor('tree.nh', 'multicolumn.fa', prior, cpd, domain)
	od=tree_order(tree)

	sequences2,tree2,factors2=generate_factor('treealt.nh', 'multicolumn.fa', prior, cpd, domain)
	od2=tree_order(tree2)
	log_likelihood=0
	log_likelihood2=0
	for i in range(1000):
		observation={}
		observation2={}
		if i%10==0:
			x=i/10+1
			stdout.write("\r%d percent" % x)
			stdout.flush()
		for k in sequences.keys():
			observation[k]=sequences[k][i]
		for l in sequences2.keys():
			observation2[l]=sequences2[l][i]
		factors_copy=copy.deepcopy(factors)
		factors_copy2=copy.deepcopy(factors2)
		log_likelihood+=VE_log(observation,tree,factors_copy,od)
		log_likelihood2+=VE_log(observation2,tree2,factors_copy2,od2)
	print
	print 'The maximal log likelihood is of tree 1 is: ',log_likelihood
	print 'The maximal log likelihood is of tree 2 is: ',log_likelihood2