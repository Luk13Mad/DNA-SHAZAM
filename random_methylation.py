from Bio import SeqIO
import random
import re
from Bio.Seq import Seq

gen_reads=SeqIO.parse("db278mtDNA.fasta","fasta")

with open("db278mtDNA_methylated.fasta", "w+") as handle:
    for read in gen_reads:
        x=random.choice(range(len(read.seq)-200))
        read.seq=read.seq[x:x+500].lower()
        
        ite=re.finditer(r"(?={})".format("c"),str(read.seq))
        pos=[i.start() for i in ite]
        idx=random.sample(pos,round(len(pos)*0.4))
        y=list(str(read.seq))
        for i in idx:
            y[i]="C"
            read.seq=Seq("".join(y))
        SeqIO.write(read, handle,"fasta")