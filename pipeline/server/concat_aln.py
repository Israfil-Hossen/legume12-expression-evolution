#!/usr/bin/env python3
# Concatenate per-orthogroup protein alignments into one supermatrix, one row per species.
# Species = header prefix before the first "_" (headers are <species>_<GeneID>).
# Usage: concat_aln.py <alignment_dir> <output.fa>
import glob, os, sys

aln_dir, out = sys.argv[1], sys.argv[2]
seqs, n, skipped = {}, 0, 0
for f in sorted(glob.glob(os.path.join(aln_dir, '*.aln'))):
    rec, name = {}, None
    for line in open(f):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split('_')[0]
            rec[name] = []
        elif name:
            rec[name].append(line)
    rec = {k: ''.join(v) for k, v in rec.items()}
    if len({len(v) for v in rec.values()}) != 1 or (seqs and set(rec) != set(seqs)):
        skipped += 1
        continue
    for k, v in rec.items():
        seqs.setdefault(k, []).append(v)
    n += 1
with open(out, 'w') as fh:
    for k in sorted(seqs):
        fh.write('>%s\n%s\n' % (k, ''.join(seqs[k])))
first = ''.join(next(iter(seqs.values())))
print('orthogroups concatenated: %d  skipped: %d  species: %d  sites: %d' % (n, skipped, len(seqs), len(first)))
