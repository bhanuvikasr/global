#!/usr/bin/env python
'''
Created on Jun 3, 2011

@author: smirarab

Edited on Jan 25, 2017

Changes made: 

TODO : convert from threshold to smallest with edges posterior
'''
import dendropy
import sys
import os
import copy
import heapq

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
if __name__ == '__main__':

    if (len(sys.argv) < 2):
        print "USAGE: %s tree_file [# edges to remove | threshold - default 75] [outfile name; - uses default] [-strip-internal|-strip-bl|strip-both|-nostrip|-remove-least-probable; default: nostrip]" %sys.argv[0]
        sys.exit(1)

    treeName = sys.argv[1]            
    t = 75 if len (sys.argv) < 3 else float(sys.argv[2])
    resultsFile="%s.%d" % (treeName,t) if len (sys.argv) < 4 or sys.argv[3]=="-" else sys.argv[3]
    
    strip_internal=True if len (sys.argv) > 4 and ( sys.argv[4]=="-strip-internal" or sys.argv[4]=="-strip-both" ) else False 
    strip_bl=True if len (sys.argv) > 4 and ( sys.argv[4]=="-strip-bl" or sys.argv[4]=="-strip-both" ) else False
    remove_least_probable=True if len (sys.argv) > 4 and ( sys.argv[4]=="-remove-least-probable") else False
    
    trees = dendropy.TreeList.get_from_path(treeName, 'newick')
    filt = lambda edge: False if (edge.label is None or (is_number(edge.label) and float(edge.label) >= (1.0*t)/100)) else True
    filt_nonleaf = lambda edge: False if (edge.label is None) else True #leaf edges have posterior as none.
    #trees_less = [trees[0]]
    for tree in trees:
        for n in tree.internal_nodes():
            if n.label is not None:
                n.label = float (n.label)
                n.edge.label = n.label
                #print n.label
                #n.label = round(n.label/2)   
        
        if remove_least_probable:
            edges_nonleaf = tree.edges(filt_nonleaf)
            least_probable_edges = heapq.nsmallest(int(t),edges_nonleaf, key=lambda e: e.label)
            print "number of least_probable_edges", len(least_probable_edges)
            for e in least_probable_edges:
                print e.label
                e.collapse()
        else:
            edges = tree.edges(filt)
            all_edges = tree.edges()
            print >>sys.stderr, len(edges), "edges will be removed", " out of", len(all_edges) 
            for e in edges:
                #print e.label
                e.collapse()
        
            if strip_internal:
                for n in tree.internal_nodes():
                    n.label = None
            if strip_bl:
                for e in tree.get_edge_set():
                    e.length = None

            #tree.reroot_at_midpoint(update_splits=False)
    print "outputting to", resultsFile
    trees.write(file=open(resultsFile,'w'),schema='newick',suppress_rooting=False)
