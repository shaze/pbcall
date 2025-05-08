#!/usr/bin/env python3

import sys
import re

f=open(sys.argv[1])
for line  in f:
    if "ngood" in line: good = int(line.split()[1])
    if "nmerr" in line: bad  = int(line.split()[1])
    if "sites_missing" in line: missing= int(line.split()[1])
total=good+bad

m = re.search(r".*trio-(\d+)-(\d+)-(\d+).rpt",sys.argv[1])
if not m:
    sys.exit(f"file name  {sys.argv[1]} not right format")

gq=m.group(1)
dp=m.group(2)
qual=m.group(3)
pct=bad*100/(good+bad)
grand_total=total+missing
g=open(sys.argv[2],"w")
g.write(f"{gq}\t{dp}\t{qual}\t{good}\t{bad}\t{total}\t{pct:.3}\t{grand_total}\n")
g.close()

