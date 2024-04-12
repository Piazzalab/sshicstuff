import os
from os.path import join, dirname
import sshicstuff as sshc

"""
1. SEQTK sub4M (subsample 4M reads) MUTANT and WT
"""

reads_dir = "/home/nicolas/Documents/Projects/sshicstuff/data/reads"
seqtk_seed = 100
seqtk_nreads = 4000000
seqtk_suffix = f"sub{seqtk_nreads // 1000000}M"


end1 = join(reads_dir, "AD162.end1.fastq.gz")
end2 = join(reads_dir, "AD162.end2.fastq.gz")
end1_out = join(reads_dir, f"AD162_{seqtk_suffix}.end1.fastq")
end2_out = join(reads_dir, f"AD162_{seqtk_suffix}.end2.fastq")

seqtk_cmd1 = f"seqtk sample -s {seqtk_seed} {end1} {seqtk_nreads} | gzip -c > {end1_out}"
seqtk_cmd2 = f"seqtk sample -s {seqtk_seed} {end2} {seqtk_nreads} | gzip -c > {end2_out}"

# apply seqtk
os.system(seqtk_cmd1)
os.system(seqtk_cmd2)

