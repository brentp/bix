"""
This compares the output of the tabix commandline tool from htslib to
this version built with the biogo bindings.
"""

import sys
import os
import subprocess
import gzip

F = sys.argv[1] if len(sys.argv) > 1 else "NA12878.wham.del.vcf.gz"

positions = []
for i, x in enumerate(gzip.open(F)):
    if x[0] == "#": continue
    if i > 2000000: break
    positions.append((x.split("\t")[0], int(x.split("\t")[1])))
print(positions[:100])
print(len(positions))

p = subprocess.Popen("go build -o main main.go", shell=True)
p.wait()

def tbxrun(query, chrom, start, end):
    p = subprocess.Popen("tabix %s %s:%d-%d" % (query, chrom, start+1, end),
                         shell=True, stdout=subprocess.PIPE)
    return p#.stdout.read().rstrip()


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
    return p, cmd#.stdout.read().rstrip(), cmd

for chrom, position in positions:
    print position
    for j in range(position-5, position+5):
        for l in range(1, 3):

            #btres = btrun(F, '1', j, j + 5)
            bxres, cmd = bxrun(F, chrom, j, j + l)
            tbres = tbxrun(F, chrom, j, j + l)
            bxres = bxres.stdout.read().strip()
            tbres = tbres.stdout.read().strip()

            #lbt = len(btres.split("\n"))
            lbx = len(bxres.split("\n")) if bxres else 0
            ltb = len(tbres.split("\n")) if tbres else 0

            assert ltb == lbx, (F, chrom, position, j, j + l, ltb, lbx, cmd)
