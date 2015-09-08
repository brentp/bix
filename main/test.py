"""
This compares the output of the tabix commandline tool from htslib to
this version built with the biogo bindings.
"""

import sys
import os
import subprocess
import gzip

F = sys.argv[1] if len(sys.argv) > 1 else "NA12878.wham.del.vcf.gz"

positions = [(x.split("\t")[0], int(x.split("\t")[1])) for x in gzip.open(F) if x[0] != "#"]
print(positions)

p = subprocess.Popen("go build -o main main.go", shell=True)
p.wait()

def tbxrun(query, chrom, start, end):
    p = subprocess.Popen("tabix %s %s:%d-%d" % (query, chrom, start+1, end),
                         shell=True, stdout=subprocess.PIPE)
    return p.stdout.read().rstrip()


def btrun(query, chrom, start, end):
    with open("t.bed", "w") as fin:
        fin.write("%s\t%d\t%d" % (chrom, start, end))
    p = subprocess.Popen("bedtools intersect -b t.bed -a %s" % query,
            shell=True, stdout=subprocess.PIPE)
    try:
        return p.stdout.read().rstrip()
    finally:
        os.unlink("t.bed")

def bxrun(query, chrom, start, end):
    cmd = "./main %s %s %d %d" % (query, chrom, start, end)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    return p.stdout.read().rstrip(), cmd

for chrom, position in positions:
    print position
    for j in range(position-10, position+10):
        for l in range(1, 10):

            #btres = btrun(F, '1', j, j + 5)
            bxres, cmd = bxrun(F, chrom, j, j + l)
            tbres = tbxrun(F, chrom, j, j + l)

            #lbt = len(btres.split("\n"))
            lbx = len(bxres.split("\n"))
            ltb = len(tbres.split("\n"))

            assert ltb == lbx, (F, chrom, position, j, j + l, ltb, lbx, cmd)
