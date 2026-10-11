#!/bin/bash
#SBATCH --job-name=is110_comprehensive
#SBATCH --output=/home/kuangh/logs/is110_comprehensive_%j.out
#SBATCH --error=/home/kuangh/logs/is110_comprehensive_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --qos=standard
#SBATCH --partition=standard

cd /home/kuangh/tools/Cycle
BASE=/groups/rubin/projects/kuang/out/IS_cycle
OUT=$BASE/is110_all_006_070
mkdir -p $OUT

echo "=== Step 1: Build consensus FASTA from Sniffles (006-070) ==="
conda run -n opfi python -c "
import csv, glob, os
csv.field_size_limit(10_000_000)
BASE = '$BASE'
OUT = '$OUT'
fasta = f'{OUT}/all_consensus.fna'
n = 0
with open(fasta, 'w') as fout:
    for b in range(6, 71):
        snif_dir = f'{BASE}/batch_{b:03d}/sniffles_output'
        if os.path.isdir(snif_dir):
            for sd in os.listdir(snif_dir):
                tf = os.path.join(snif_dir, sd, f'{sd}.table.txt')
                if not os.path.exists(tf): continue
                with open(tf) as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        seq = row.get('Consensus', '')
                        if len(seq) >= 500:
                            fout.write(f'>{row[\"UUID\"]}\n{seq}\n')
                            n += 1
        else:
            for cdir in sorted(glob.glob(f'{BASE}/batch_{b:03d}/chunks/chunk_*/sniffles_output')):
                for sd in os.listdir(cdir):
                    tf = os.path.join(cdir, sd, f'{sd}.table.txt')
                    if not os.path.exists(tf): continue
                    with open(tf) as f:
                        reader = csv.DictReader(f, delimiter='\t')
                        for row in reader:
                            seq = row.get('Consensus', '')
                            if len(seq) >= 500:
                                fout.write(f'>{row[\"UUID\"]}\n{seq}\n')
                                n += 1
print(f'Total sequences: {n}')
"

echo "=== Step 2: Parallel Prodigal (48 chunks) ==="
mkdir -p $OUT/prodigal_chunks
conda run -n opfi python -c "
fasta = '$OUT/all_consensus.fna'
out_dir = '$OUT/prodigal_chunks'
seqs = []
with open(fasta) as f:
    header = None; seq_lines = []
    for line in f:
        if line.startswith('>'):
            if header: seqs.append((header, ''.join(seq_lines)))
            header = line; seq_lines = []
        else: seq_lines.append(line)
    if header: seqs.append((header, ''.join(seq_lines)))
print(f'Total: {len(seqs)} sequences')
chunk_size = (len(seqs) + 48 - 1) // 48
for i in range(48):
    chunk = seqs[i*chunk_size:(i+1)*chunk_size]
    if not chunk: continue
    with open(f'{out_dir}/chunk_{i}.fna', 'w') as f:
        for h, s in chunk: f.write(h + s)
"

for i in $(seq 0 47); do
    CHUNK=$OUT/prodigal_chunks/chunk_${i}.fna
    [ ! -f "$CHUNK" ] && continue
    conda run -n opfi prodigal -i $CHUNK -a $OUT/prodigal_chunks/chunk_${i}.faa -o /dev/null -p meta -q &
done
wait
cat $OUT/prodigal_chunks/chunk_*.faa > $OUT/all_proteins.faa
echo "Proteins: $(grep -c '^>' $OUT/all_proteins.faa)"

echo "=== Step 3: HMM ==="
conda run -n opfi python -c "
import json, logging, os, sys
from collections import defaultdict
sys.path.insert(0, '/home/kuangh/tools/Cycle')
from cycle.is110_filter import IS110Filter
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')

OUT = '$OUT'
uuid_orfs = defaultdict(list)
with open(f'{OUT}/all_proteins.faa') as f:
    header = None; seq_lines = []
    for line in f:
        if line.startswith('>'):
            if header and seq_lines:
                seq = ''.join(seq_lines).rstrip('*')
                uuid_orfs[header['uuid']].append({'start': header['start'], 'end': header['end'], 'strand': header['strand'], 'sequence': seq})
            parts = line[1:].strip().split(' # ')
            if len(parts) < 4:
                header = None; seq_lines = []; continue
            prot_id = parts[0]
            uuid = '_'.join(prot_id.rsplit('_', 1)[:-1])
            strand = '+' if parts[3] == '1' else '-'
            header = {'uuid': uuid, 'start': int(parts[1]), 'end': int(parts[2]), 'strand': strand}
            seq_lines = []
        else: seq_lines.append(line.strip())

hmm_fasta = f'{OUT}/all_proteins_hmm.faa'
with open(hmm_fasta, 'w') as f:
    for uuid, orfs in uuid_orfs.items():
        for orf in orfs:
            f.write(f'>{uuid}__{orf[\"start\"]}_{orf[\"end\"]}_{orf[\"strand\"]}\n{orf[\"sequence\"]}\n')

filt = IS110Filter()
is110_ids, _ = filt.run_hmmsearch(hmm_fasta, OUT, cpus=44)
with open(f'{OUT}/is110_ids.json', 'w') as f:
    json.dump(list(is110_ids), f)
print(f'IS110 IDs: {len(is110_ids)}')
"

echo "=== Done ==="
