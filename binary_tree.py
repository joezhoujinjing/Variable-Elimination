###############################################################################
# simple binary tree class.  Will be used to construct the clique tree.
#
# author: Jinjing Zhou, Isaac Caswell
# date: Jan 15th, 2016
###############################################################################

class BinaryTree:
    def __init__(self):
        self.idx = None # self.idx is an key/index
        self.data = None # self.data will be name or data
        self.left = None # self.left is a BinaryTree object
        self.right = None # self.right is also a BinaryTree object.


    def is_leaf(self):
        return self.left == None and self.right == None

    def parse_from_string(self, str_to_parse):
        """
        parses a BinaryTree from a string of nested parentheses.    
        The string str_to_parse might look like the following:

        ((((((human,baboon),marmoset),((rat,mouse),rabbit)),(snail,goat)),elephant),platypus)
        """

        if '(' not in str_to_parse and ',' not in str_to_parse:
            # it's a leaf node! 
            self.data = str_to_parse
            return

        str_to_parse = str_to_parse[1:-1] #remove parentheses
        num_open = 0
        for i in xrange(len(str_to_parse)):
            if num_open == 0 and str_to_parse[i] == ',':
                self.left = BinaryTree()
                self.left.parse_from_string(str_to_parse[:i])
                self.right = BinaryTree()
                self.right.parse_from_string(str_to_parse[(i+1):])
                break
            elif str_to_parse[i] == '(':
                num_open += 1
            elif str_to_parse[i] == ')':
                num_open -= 1

    def __repr__(self):
        """
        Returns a representation of the binary tree.    This may be useful for debugging.
        usage:

        tree = BinaryTree
        print repr(tree)

        A note for the curious: obj.__repr__(), corresponding to the builtin function repr(obj),
        and obj.__str__(), corresponding to the builtin function str(obj), are identical except 
        that the former aims for being informative (and therefore often longer) whereas the latter 
        aims for being readable.
        """
        if self.left == None and self.right == None: # a leaf node
            return self.data
        return '(' + repr(self.left) + ',' + repr(self.right) + ')'

    def assign_idx(self, interiorNodes, rf):
        """
        :param list rf: a singleton list to keep track of the indices. It should by all rights be 
            an integer, but it needs to be referenced by all the various stack frames, so this was 
            an easy workaround.

        This function labels each interior node with an index corresponding to its place in the tree.
        For instance, in a tree with three leaf nodes and two interior nodes, the interior nodes 
        will be labeled "animal_1" and "animal_0"

        """
        if self.left != None:
            self.left.assign_idx(interiorNodes, rf)
        if self.right != None:
            self.right.assign_idx(interiorNodes, rf)
        self.idx = rf[0]
        if self.data == None:
            interiorNodes.append(self.idx)
            self.data = 'animal_' + str(self.idx+1)
        rf[0] += 1
