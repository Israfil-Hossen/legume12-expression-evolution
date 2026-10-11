#!/usr/bin/env python3
# Step 3A: reciprocal-best-hit paralog pairs, Ks with MAFFT + pal2nal + yn00, MCScanX class.
# Usage: pairs_ks.py <species> <threads>
# Needs on PATH: mafft, pal2nal.pl, yn00. Reads step3/<sp>/{<sp>.blast, proteins.fa, cds.fa,
# <sp>.collinearity, <sp>.gene_type}. Writes step3/<sp>/pairs_ks.tsv.
import os, subprocess, sys, shutil
from concurrent.futures import ThreadPoolExecutor

P = os.environ.get('LEGUME2', '/home/liuchuanfa/israfil/legume2')
sp, threads = sys.argv[1], int(sys.argv[2])
d = f'{P}/step3/{sp}'

def read_fa(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}

prot, cds = read_fa(f'{d}/proteins.fa'), read_fa(f'{d}/cds.fa')

best = {}                                   # query -> (bitscore, subject), self hits ignored
for line in open(f'{d}/{sp}.blast'):
    f = line.rstrip('\n').split('\t')
    q, s, bits = f[0], f[1], float(f[11])
    if q != s and (q not in best or bits > best[q][0]):
        best[q] = (bits, s)
pairs = sorted({tuple(sorted((q, s))) for q, (b, s) in best.items() if best.get(s, (0, None))[1] == q})

collinear = set()
for line in open(f'{d}/{sp}.collinearity'):
    if line.startswith('#') or not line.strip():
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) >= 3:
        collinear.add(tuple(sorted((f[1].strip(), f[2].strip()))))
gtype = {}
for line in open(f'{d}/{sp}.gene_type'):
    f = line.split()
    if len(f) == 2:
        gtype[f[0]] = f[1]                  # 0 singleton, 1 dispersed, 2 proximal, 3 tandem, 4 WGD/segmental

work = f'{d}/ks_work'
os.makedirs(work, exist_ok=True)
CTL = 'seqfile = codon.phy\noutfile = yn\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n'

def ks(i_pair):
    i, (a, b) = i_pair
    w = f'{work}/p{i}'
    try:
        if a not in prot or b not in prot or a not in cds or b not in cds:
            return None
        os.makedirs(w, exist_ok=True)
        open(f'{w}/prot.fa', 'w').write(f'>{a}\n{prot[a]}\n>{b}\n{prot[b]}\n')
        open(f'{w}/cds.fa', 'w').write(f'>{a}\n{cds[a]}\n>{b}\n{cds[b]}\n')
        aln = subprocess.run(['mafft', '--quiet', '--auto', 'prot.fa'], cwd=w, capture_output=True, text=True).stdout
        open(f'{w}/prot.aln', 'w').write(aln)
        phy = subprocess.run(['pal2nal.pl', 'prot.aln', 'cds.fa', '-output', 'paml', '-nogap'],
                             cwd=w, capture_output=True, text=True).stdout
        if not phy.strip():
            return None
        open(f'{w}/codon.phy', 'w').write(phy)
        open(f'{w}/yn00.ctl', 'w').write(CTL)
        subprocess.run(['yn00', 'yn00.ctl'], cwd=w, capture_output=True, text=True)
        last = [l for l in open(f'{w}/2YN.dS') if l.strip()][-1].split()
        return (a, b, float(last[-1]))
    except Exception:
        return None
    finally:
        shutil.rmtree(w, ignore_errors=True)

with ThreadPoolExecutor(max_workers=threads) as ex:
    res = [r for r in ex.map(ks, enumerate(pairs)) if r]

with open(f'{d}/pairs_ks.tsv', 'w') as fo:
    fo.write('gene1\tgene2\tks\tcollinear\ttype1\ttype2\n')
    for a, b, k in res:
        fo.write(f'{a}\t{b}\t{k}\t{int((a, b) in collinear)}\t{gtype.get(a, "NA")}\t{gtype.get(b, "NA")}\n')
print(f'{sp}: reciprocal-best-hit pairs {len(pairs)}; with Ks {len(res)}; collinear {sum((a, b) in collinear for a, b, _ in res)}')
