#!/usr/bin/env python
'''
Created on Jun 3, 2011

@author: smirarab

Edited on Jan 25, 2017

Changes made: added new feature to remove least probable branches: access with -remove-least-probable
              and now threshold is divided by 100 before using. i.e.  input 75 now means 0.75

TODO : take threshold values in a file and use it.
'''
import dendropy
import sys
import os
import copy
import heapq
import csv
# import pandas

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

    trees = dendropy.TreeList.get_from_path(treeName, 'newick')

    if len (sys.argv) < 3:
        t_list = [75]*len(trees)
    elif is_number(sys.argv[2]):
        t_list  = [float(sys.argv[2])]*len(trees)
    else: # read from the file. If it fails then break.
        # print "reading thresholds from the file"
        with open(sys.argv[2], 'rb') as f:
            reader = csv.reader(f)
            l = list(reader)
            t_list = [float(x[0]) if is_number(x[0]) else 0.0 for x in l]
        print "Finished reading. Made a list of length", len(t_list) #, t_list[553674], t_list[553673], t_list[553672]
        if len(t_list)<len(trees):
            print "Not enough threshold values"

    resultsFile="%s.%d" % (treeName,t) if len (sys.argv) < 4 or sys.argv[3]=="-" else sys.argv[3]

    strip_internal=True if len (sys.argv) > 4 and ( sys.argv[4]=="-strip-internal" or sys.argv[4]=="-strip-both" ) else False
    strip_bl=True if len (sys.argv) > 4 and ( sys.argv[4]=="-strip-bl" or sys.argv[4]=="-strip-both" ) else False
    remove_least_probable=True if len (sys.argv) > 4 and ( sys.argv[4]=="-remove-least-probable") else False

    if remove_least_probable:
        filt_nonleaf = lambda edge: False if (edge.label is None) else True #leaf edges have posterior as none.
    else:
        filt = lambda edge: False if (edge.label is None or (is_number(edge.label) and float(edge.label) >= (1.0*t)/100)) else True #try if this works with different t values.

    #trees_less = [trees[0], trees[1]]
    #err = [t_list[0],t_list[1]]
    for (tree, t) in zip(trees,t_list):
    #for tree in trees:
        for n in tree.internal_nodes():
            if n.label is not None:
                n.label = float (n.label)
                n.edge.label = n.label
                #print n.label
                #n.label = round(n.label/2)

        if remove_least_probable:
            edges_nonleaf = tree.edges(filt_nonleaf)
            # print t, len(edges_nonleaf), max(int(2.0*t*len(edges_nonleaf)),len(edges_nonleaf)
            least_probable_edges = heapq.nsmallest(min(int(round(1.0*t*len(edges_nonleaf))),len(edges_nonleaf)),edges_nonleaf, key=lambda e: e.label)
            #print len(edges_nonleaf)

            # print "number of least_probable_edges", len(least_probable_edges)
            for e in least_probable_edges:
                # print e.label
                e.collapse()
        else:
            edges = tree.edges(filt)
            all_edges = tree.edges()
            # print >>sys.stderr, len(edges), "edges will be removed", " out of", len(all_edges)
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
