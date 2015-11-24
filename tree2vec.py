import dendropy
import scipy.sparse

def tree2vec(tree, tn=None):
    if isinstance(tree, dendropy.Tree):
        t = tree
        tn = t.taxon_namespace
    elif tree[0] == '(':
        t = dendropy.Tree.get_from_string(tree, 'newick', taxon_namespace = tn)
    else:
        t = dendropy.Tree.get_from_path(tree, 'newick', taxon_namespace = tn)

    t.encode_bipartitions()

    outv = []
    
    for b in t.bipartition_encoding:
        outv.append(b.leafset_bitmask)
        
    return outv 
