#! /usr/bin/env python
 
import os
import sys
import re
import subprocess as sub

name_pt = re.compile("(?<=[>])(.+)")
 
if ("--help" in sys.argv) or ("-?" in sys.argv) or len(sys.argv) < 1:
    sys.stderr.write("usage: %s <alignment-file-path> [<results-file-path>]\n"%sys.argv[0])
    sys.exit(1)
 
src_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
if not os.path.exists(src_fpath):
    sys.stderr.write('Not found: "%s"' % src_fpath)

p = sub.Popen(['simplifyfasta.sh',src_fpath],stdout=sub.PIPE,stderr=sub.PIPE)
simpl, errors = p.communicate()
if errors is not None and errors != "":
    print errors
    sys.exit(1)

dest_fpath=None
if len(sys.argv) > 2:
    dest_fpath=os.path.expanduser(os.path.expandvars(sys.argv[2]))
    dest = open(dest_fpath, "w")
else:
    dest = sys.stdout

print >>dest, "%s %s %s %s %s %s %s %s %s %s" %("SEQUENCE","TAXON","A_C","C_C","G_C","T_C","A_R","C_R","G_R","T_R")        
#print mapping
try:
    seq=""
    for l in simpl.split("\n"):    
        if l.startswith(">"):
            seq=l[1:]
            seq2=seq.split("_")[0]
        else:
            a= [0, 0, 0, 0]
            c= [0, 0, 0, 0]
            g= [0, 0, 0, 0]
            t= [0, 0, 0, 0]
            for x in l: 
                a[0]+= (1 if x == "A" or x == "a" else 0)  
                c[0]+= (1 if x == "C" or x == "c" else 0)
                g[0]+= (1 if x == "G" or x == "g" else 0)
                t[0]+= (1 if x == "T" or x == "t" else 0)
            
            s=a[0]+c[0]+g[0]+t[0]+0.0
            sc=s/3
            
            lc = l[0::3]
            for x in lc: 
                a[1]+= (1 if x == "A" or x == "a" else 0)  
                c[1]+= (1 if x == "C" or x == "c" else 0)
                g[1]+= (1 if x == "G" or x == "g" else 0)
                t[1]+= (1 if x == "T" or x == "t" else 0)
           
            lc = l[1::3]
            for x in lc: 
                a[2]+= (1 if x == "A" or x == "a" else 0)  
                c[2]+= (1 if x == "C" or x == "c" else 0)
                g[2]+= (1 if x == "G" or x == "g" else 0)
                t[2]+= (1 if x == "T" or x == "t" else 0)
            
            lc = l[2::3]
            for x in lc: 
                a[3]+= (1 if x == "A" or x == "a" else 0)  
                c[3]+= (1 if x == "C" or x == "c" else 0)
                g[3]+= (1 if x == "G" or x == "g" else 0)
                t[3]+= (1 if x == "T" or x == "t" else 0)
            
            print >>dest, "%s %s %d %d %d %d %f %f %f %f %d %d %d %d %f %f %f %f %d %d %d %d %f %f %f %f %d %d %d %d %f %f %f %f" %(
                          seq,seq2,
                                  a[0],c[0],g[0],t[0],a[0]/s,c[0]/s,g[0]/s,t[0]/s,
                                                         a[1],c[1],g[1],t[1],a[1]/sc,c[1]/sc,g[1]/sc,t[1]/sc,
                                                                                  a[2],c[2],g[2],t[2],a[2]/sc,c[2]/sc,g[2]/sc,t[2]/sc,
                                                                                                        a[3],c[3],g[3],t[3],a[3]/sc,c[3]/sc,g[3]/sc,t[3]/sc)
except RuntimeError as e:
    dest.close()
    if dest_fpath is not None:
        os.remove(dest_fpath)
    raise e
        
dest.close()
