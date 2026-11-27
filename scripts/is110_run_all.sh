#!/bin/bash
#SBATCH --job-name=is110_all
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G

# Site settings (partition/qos/account, log directory, data and asset paths) come from config/.
# Submit with:  bin/cycle-sbatch scripts/is110_run_all.sh
source "${CYCLE_ROOT:?submit with bin/cycle-sbatch or export CYCLE_ROOT}/config/load_site.sh"

cd "$CYCLE_ROOT"
BASE=$CYCLE_DATA_ROOT
OUT=$BASE/is110_all_006_052
mkdir -p $OUT

echo "=== Step 1: Build consensus FASTA from Sniffles ==="
conda run -p "$CYCLE_ENV_PREFIX" python -c "
import csv, glob, os
csv.field_size_limit(10_000_000)
BASE = '$BASE'
OUT = '$OUT'
fasta = f'{OUT}/all_consensus.fna'
n = 0
with open(fasta, 'w') as fout:
    for b in range(6, 53):
        found = False
        # Try top-level sniffles_output
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
            found = True
        # Try chunk-level sniffles_output
        if not found:
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
# Chunks of an earlier run must not be merged into this one.
rm -rf "$OUT/prodigal_chunks"
mkdir -p $OUT/prodigal_chunks
conda run -p "$CYCLE_ENV_PREFIX" python -c "
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

pids=""
for i in $(seq 0 47); do
    CHUNK=$OUT/prodigal_chunks/chunk_${i}.fna
    [ ! -f "$CHUNK" ] && continue
    conda run -p "$CYCLE_ENV_PREFIX" prodigal -i $CHUNK -a $OUT/prodigal_chunks/chunk_${i}.faa -o /dev/null -p meta -q &
    pids="$pids $!"
done
for pid in $pids; do
    wait $pid || { echo "A Prodigal worker failed" >&2; exit 1; }
done
cat $OUT/prodigal_chunks/chunk_*.faa > $OUT/all_proteins.faa
echo "Proteins: $(grep -c '^>' $OUT/all_proteins.faa)"

echo "=== Step 3: Parse + HMM ==="
conda run -p "$CYCLE_ENV_PREFIX" python -c "
import json, logging, os, sys
from collections import defaultdict
sys.path.insert(0, '$CYCLE_ROOT')
from cycle.is110_filter import IS110Filter
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

OUT = '$OUT'
uuid_orfs = defaultdict(list)
n_proteins = 0
with open(f'{OUT}/all_proteins.faa') as f:
    header = None; seq_lines = []
    for line in f:
        if line.startswith('>'):
            if header and seq_lines:
                seq = ''.join(seq_lines).rstrip('*')
                uuid_orfs[header['uuid']].append({'start': header['start'], 'end': header['end'], 'strand': header['strand'], 'sequence': seq})
                n_proteins += 1
            parts = line[1:].strip().split(' # ')
            if len(parts) < 4:
                header = None; seq_lines = []; continue
            prot_id = parts[0]
            uuid = '_'.join(prot_id.rsplit('_', 1)[:-1])
            strand = '+' if parts[3] == '1' else '-'
            header = {'uuid': uuid, 'start': int(parts[1]), 'end': int(parts[2]), 'strand': strand}
            seq_lines = []
        else: seq_lines.append(line.strip())
    if header and seq_lines:
        seq = ''.join(seq_lines).rstrip('*')
        uuid_orfs[header['uuid']].append({'start': header['start'], 'end': header['end'], 'strand': header['strand'], 'sequence': seq})
        n_proteins += 1
logger.info('Parsed %d proteins from %d elements', n_proteins, len(uuid_orfs))

hmm_fasta = f'{OUT}/all_proteins_hmm.faa'
with open(hmm_fasta, 'w') as f:
    for uuid, orfs in uuid_orfs.items():
        for orf in orfs:
            f.write(f'>{uuid}__{orf[\"start\"]}_{orf[\"end\"]}_{orf[\"strand\"]}\n{orf[\"sequence\"]}\n')

filt = IS110Filter()
is110_ids, trans_domains = filt.run_hmmsearch(hmm_fasta, OUT, cpus=44)
with open(f'{OUT}/is110_ids.json', 'w') as f:
    json.dump(list(is110_ids), f)
logger.info('Saved %d IS110 IDs', len(is110_ids))
"

echo "=== Done ==="
