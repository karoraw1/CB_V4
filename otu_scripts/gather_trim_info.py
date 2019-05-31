import os 
# run this from inside processing dir
dirs = []
for i in os.listdir(os.getcwd()):
    if os.path.exists(os.path.join(i, 'FASTQ')):
        dirs.append(i)

all_fds = []
for d in dirs:
     dd = os.path.join(d, 'FASTQ')
     ddd = os.path.join(dd, "Demux")
     ddt = os.path.join(dd, 'Trim')
     all_fds += [(i[:-9], os.path.join(ddd, i), os.path.join(ddt, i[:-9]+"_F_filt.fastq")) for i in os.listdir(ddd) if i.endswith('R1.fastq')]

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def count_lines2(filepath):
    with open(filepath, "r", encoding="utf-8", errors='ignore') as f:
        count = sum(bl.count("\n") for bl in blocks(f))
    return count/4

def all_counts(arg):
    x, y, z = arg
    if os.path.exists(y):
        rc = count_lines2(y)
    else:
        rc = 'NA'
    if os.path.exists(z):
        tc = count_lines2(z)
    else:
        tc = 'NA'
    print(x, rc, tc)
    return (x, rc, tc)

existers = map(all_counts, all_fds)

to_write = "\n".join(["\t".join([i, str(j), str(k)]) for i, j, k in existers])

f_out = "read_counts.tsv"

with open(f_out, 'w') as foh:
    foh.write(to_write)
