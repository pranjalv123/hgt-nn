import dendropy
import subprocess
import numpy
from ASTRID import ASTRID

#you need ASTRID, which can be installed with
# pip install astrid-phylo
#and fasttree, which can be gotten from apt-get or at
#http://meta.microbesonline.org/fasttree/

def get_desc_edges(t, n):
    desc = set()
    nodes = [n]
    while nodes:
        n = nodes.pop()
        desc.add(n.edge)
        for n2 in n.child_nodes():
            nodes.append(n2)
    return list(desc)

class Simhgt(object):
    def __init__(self, stree):
        self.stree = stree
        self.stree.suppress_unifurcations()
        
    def dohgt(self, ltree, e1, e2):

        #prune
        length = e1.length
        h = e1.head_node
        ltree.prune_subtree(h, suppress_unifurcations=False)
        h.edge_length=length

        #regraft
        tn = e2.tail_node
        hn = e2.head_node
        new = tn.new_child()
        new.edge.length = e2.length - length
        tn.remove_child(hn)
        new.add_child(hn)
        hn.edge.length = length
        new.add_child(h)
        ltree.update_splits(suppress_unifurcations=True)
        return ltree

    def dorandomhgt(self):
        ltree = self.stree.clone()        
        
        top  = numpy.random.choice(ltree.internal_nodes(), 1)[0]

        ltree.update_bipartitions()
        
        e1 = numpy.random.choice(get_desc_edges(ltree, top.child_nodes()[0]), 1)[0]
        e2 = numpy.random.choice(get_desc_edges(ltree, top.child_nodes()[1]), 1)[0]

        c1 = e1.split_bitmask
        c2 = e2.split_bitmask
        
        ltree = self.dohgt(ltree, e1, e2)
        for edge in ltree.postorder_edge_iter():
            edge.pop_size = 100000
        return ltree, c1, c2

    def simgtrees(self,  seqlen=100, nhgt=10, nils=10):

        gtrees = dendropy.TreeList()
        hgtevents = []
        for i in range(nhgt):
            ltree, c1, c2  = self.dorandomhgt()
            gtrees.append(ltree)
            
            for i in range(nils):
                gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(containing_taxon_namespace=ltree.taxon_namespace, num_contained=1, contained_taxon_label_fn = lambda a,b:a)
                gtree = dendropy.simulate.treesim.contained_coalescent_tree(ltree, gene_to_containing_taxon_map=gene_to_species_map)
                gtrees.append(gtree)
                hgtevents.append((c1, c2))

        seqs = []

        for tree in gtrees:
            seqs.append(dendropy.simulate.charsim.hky85_chars(seqlen, tree))

        self.truetrees = gtrees
        self.seqs = seqs
        self.hgtevents = hgtevents
        return gtrees, seqs, hgtevents

    def estgtrees(self):
        self.esttrees = dendropy.TreeList()
        for mat in self.seqs:
            tree, _ = subprocess.Popen(["fasttree", "-nt", "-gtr"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate(mat.as_string(schema='phylip'))
            self.esttrees.append(dendropy.Tree.get_from_string(tree, 'newick'))
            
            
def gen_test_data(file=None, seed=None, nreps=50, ntaxa=15, ngenes=1000):
    if seed:
        numpy.random.seed(seed)
        dendropy.utility.GLOBAL_RNG.seed(seed)
    if file:
        sh = Simhgt(dendropy.Tree.get_from_file(file, 'newick'))
    else:
        stree = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=ntaxa)
        sh = Simhgt(stree)


    egts = []
    tgts = []
    seqs = []
    hgts = []
    
    for i in range(nreps):
        sh.simgtrees(nhgt=ngenes, nils=1)
        sh.estgtrees()
        egts.append(sh.esttrees)
        tgts.append(sh.truetrees)
        seqs.append(sh.seqs)
        hgts.append(sh.hgtevents)
    return tgts, egts, seqs, hgts, sh.stree



def pipeline(truetrees, esttrees,  seqs, hgtevents, stree, learner=None):
    ast = ASTRID(esttrees)
    ast.run('auto')
    estspeciestree = ast.tree

    print str(stree)
    
    estspeciestree.migrate_taxon_namespace(stree.taxon_namespace)

    print str(estspeciestree)
    
    print "RF distance", (dendropy.treecalc.treecompare.false_positives_and_negatives(estspeciestree, stree)), len(stree.internal_edges())

    sim = Simhgt(estspeciestree)
    
    
    sim.simgtrees(seqlen=20, nhgt=10, nils=10)
    sim.estgtrees()
    if not learner:
        return
    traintrees = sim.truetrees + sim.esttrees
    trainhgts = sim.hgtevents + sim.hgtevents
    
    learner.train(traintrees, trainhgts)

    estimated_hgts = learner.get(esttrees)
    return estimated_hgts


    
