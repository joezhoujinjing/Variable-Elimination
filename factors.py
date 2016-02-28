###############################################################################
# utility functions for manipulating factors; ported from Daphne Koller's 
# Matlab utility code
# author: Jinjing Zhou, Billy Jun
# date: Jan 20, 2016
###############################################################################

import numpy as np
import copy

def assignment_to_indices(A, card):
    """
    :param - A: an assignment
    :param list card: a list of the cardinalities of the variables in the assignment
    """
    A = np.array(A, copy=False)
    card = np.array(card, copy=False)
    C = card.flatten()
    if np.any(np.shape(A) == 1):
        I = np.cumprod(np.concatenate(([1.0], C[:0:-1]))) * (A.T).flatten()
    else:
        B = A[:,::-1]
        I = np.sum(np.tile(np.cumprod(np.concatenate(([1.0], C[:0:-1]))), \
                (B.shape[0], 1)) * B, axis=1)
    return np.array(I, dtype='int32')

def indices_to_assignment(I, card):
    """
    :param - I: a list of indices
    :param list card: a list of the cardinalities of the variables in the assignment
    """
    I = np.array(I, copy=False)
    card = np.array(card, copy=False)
    C = card.flatten()
    A = np.mod(np.floor(
            np.tile(I.flatten().T, (len(card), 1)).T / \
            np.tile(np.cumprod(np.concatenate(([1.0], C[:0:-1]))), (len(I), 1))), \
            np.tile(C[::-1], (len(I), 1)))
    return A[:,::-1]


def intersection_indices(a, b):
    """
    :param list a, b: two lists of variables from different factors.

    returns a tuple of 
        (indices in a of the variables are in the intersection of a.scope and b.scope
indices of those same variables within the list b)
    """
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    mapA = []
    mapB = []
    for i, itm in enumerate(a):
        if itm in bind:
            mapA.append(i)
            mapB.append(bind.get(itm))
    return mapA, mapB

class Factor:
    def __init__(self, f=None, scope=[], card=[], val=None, name=None):
        """
        :param Factor f: if this parameter is not None, then the constructor makes a 
            copy of it.
        :param list scope: a list of variable names that are in the scope of this factor
        :param list card: a list of integers coresponding to the cardinality of each variable
            in scope
        :param np.ndarray val: an array coresponding to the values of different assignments 
            to the factor. val is a numpy.ndarray of shape self.card. Therefore, if this factor is over
            three binary variables, self.val will be an array of shape (2,2,2)
        :param str name: the name of the factor.  Useful for debugging only--no functional
            purpose.
        """
        assert len(scope) == len(card)

        # self.scope: a list of the variables over which this Factor defines a ddistribution
        self.scope = scope

        # self.card: the cardinality of each variable in self.scope
        self.card = card

        # use the name field for debugging, imo
        self.name = name

        self.val = val

        if f is not None:
            self.scope = list(f.scope)
            self.card = list(f.card)
            self.val = np.array(f.val, copy=True)
            self.name = f.name

    def compose_factors(self, f, operator, opname="op"):

        g = Factor() # modify this to be the composition of two Factors and then return it
        #g.name = "(%s %s %s)"%(self.name, opname, f.name)
        ###############################################################################
        #add self and f
        new_scope=list(set(copy.deepcopy(self.scope)+copy.deepcopy(f.scope)))
        new_card=[]
        for s in new_scope:
            if s in self.scope:
                new_card.append(self.card[self.scope.index(s)])
            else:
                new_card.append(f.card[f.scope.index(s)])
        g.scope=new_scope
        g.card=new_card
        g.val=np.ones(g.card)


        assign=indices_to_assignment(range(np.prod(g.card)),g.card)
        for a in assign:         
            a_s=()
            a_f=()
            for s_s in self.scope:
                a_s+=(a[g.scope.index(s_s)],)
            for f_s in f.scope: 
                a_f+=(a[g.scope.index(f_s)],)
            g.val[tuple(a)]=operator(self.val[a_s],f.val[a_f])
        return g

    def sum(self, f):
        """
        Returns a factor that is the result of adding this factor with factor f.

        :param Factor f: the factor by which to multiply this factor.
        :rtype: Factor
        """
        return self.compose_factors(f, operator=lambda x, y: x+y, opname = "+")   


    def multiply(self, f):
        """
        Returns a factor that is the result of multiplying this factor with factor f.
        :param Factor f: the factor by which to multiply this factor.
        :rtype: Factor
        """
        pass
        ###############################################################################
        # TODO: Your code here! 
        return self.compose_factors(f, operator=lambda x, y: x*y, opname = "*")              

    def marginalize_all_but(self, var, marginal_type="sum"):
        """
        returns a factor that is like unto this one except that all variables except those 
        in the set var have been marginalized out.

        :param set var: a set of the variables not to be marginalized out.
        :param str marginal_type: either "sum", signifying sum-marginalization,
            or  "max", signifying max-marginalization.
        :rtype: Factor 

        """
        assert marginal_type=="sum" or marginal_type=="max"
        g = Factor() 
        marginalized_out = ", ".join([str(v) for v in set(self.scope) - set(var)])
        #g.name = "(\sum_{%s} %s)"%(marginalized_out, self.name)
        g.scope=list(set(var))
        g.card=[]
        for s in g.scope:
            g.card.append(self.card[self.scope.index(s)])
        g.val=np.zeros(g.card)

        assign=indices_to_assignment(range(np.prod(self.card)),self.card)
        
        if marginal_type=='sum':
            for a in assign:
                g_s=()
                for s in g.scope:
                    g_s+=(a[self.scope.index(s)],)
                g.val[g_s]+=self.val[tuple(a)]
        
        if marginal_type=='max':
            for a in assign:
                g_s=()
                for s in g.scope:
                    g_s+=(a[self.scope.index(s)],)
                g.val[g_s]=max(g.val[g_s],self.val[tuple(a)])
        
        return g   
  

    def observe(self, var, val):
        """
        Returns a version of this factor with variable var observed as having taken on value val.
        if var is not in the scope of this Factor, a duplicate of this factor is returned.
        
        :param str var: the observed variable
        :param int val: the value that variable took on
        :return: a Factor corresponding to this factor with var observed at val
        """
        f = Factor(self) #make copy.  You'll modify this.
        #f.name = "(%s with variable %s observed as %s)"%(self.name, var, val)  
        assign=indices_to_assignment(range(np.prod(f.card)),f.card)
        i = f.scope.index(var)
        for a in assign:
            if a[i]!=val:
                f.val[tuple(a)]=0
        return f

    def __repr__(self):
        """
        returns a descriptive string representing this factor!
        """
        r = "Factor object with scope %s and corresponding cardinalities %s"%(self.scope, self.card)
        r += "\nCPD:\n" + str(self.val)
        if self.name:
            r = "Factor %s:\n"%self.name + r
        return r + "\n"

    def __str__(self):
        """
        returns a nice string representing this factor!  Note that we can now use string formatting
        with %s and this will cast our class into somethign nice and readable.
        """
        return self.name   



