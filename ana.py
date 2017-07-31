#!/usr/bin/env python

import sys
file_name = sys.argv[1]
with open(file_name,"r") as f:
    data = f.read().split("\n")[-3]
    f.close()

L = map(int,file_name.split("."))
P = L[0]*L[1]

dd = [(bin(i)[2:].zfill(P),j) for i,j in enumerate(map(float,data.split(" ")[:-1])) if j != 0]

dd.sort(key=lambda x:x[1])

def f(d):
    ans = ""
    for i in xrange(L[0]):
        ans += d[:L[1]]+"\n"
        d = d[L[1]:]
    return ans

with open("%s.data"%file_name,"w") as ff:
    for i in dd:
        print >>ff, f(i[0]),i[1],"\n"
