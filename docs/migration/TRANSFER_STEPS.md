# Moving the verified bundle from IGI to Lawrencium, step by step

Everything lands in one folder on Lawrencium, here `TARGET=$HOME/IS110_ecology/Cycle_from_IGI`:

```
$TARGET/
  transfer/   the archives exactly as copied from IGI (~438 GB)
  data/       unpacked: IS_cycle/, cycle_rna_analysis/, IS110_downstream/, ECOR_batch/, assets/ (~1.28 TB)
  code/Cycle  this repository, branch migration/lawrencium
```

Peak space is about 1.8 TB. The IS_cycle reference genomes are restored as hard links, so `data/`
must be on a filesystem that supports them (Lawrencium scratch does). Replace `<user>`, `<account>`,
`<partition>` and `<qos>` with your own values; nothing private belongs in this file.

## 1. Lawrencium: check space, tools and Slurm account

```bash
realpath ~/IS110_ecology
df -h "$(realpath ~/IS110_ecology)"
sacctmgr show assoc user=$USER format=account%30,partition,qos%40
zstd --version; python3 --version     # zstd >= 1.4 and python3 >= 3.6 are enough (archives use 2 MiB windows)
mkdir -p ~/logs
```

If `~/IS110_ecology` is on the home filesystem (30 GB quota), keep the folder on scratch and link it
in, so it still appears as `~/IS110_ecology/Cycle_from_IGI`:

```bash
mkdir -p /global/scratch/users/$USER/Cycle_from_IGI/transfer
ln -s /global/scratch/users/$USER/Cycle_from_IGI ~/IS110_ecology/Cycle_from_IGI
```

Otherwise (enough free space where it is): `mkdir -p ~/IS110_ecology/Cycle_from_IGI/transfer`.
Scratch is not backed up and is subject to a purge policy, so keep the IGI sources until the
archives in `transfer/` also have a permanent home.

## 2. IGI: transfer (run it yourself; asks for the Lawrencium PIN + one-time password)

```bash
tmux new -s xfer
cd /groups/rubin/projects/kuang/migration_20261117
export DEST=<user>@lrc-xfer.lbl.gov:/global/home/users/<user>/IS110_ecology/Cycle_from_IGI/transfer
bash migrate.sh --dry-run      # lists what would be copied
bash migrate.sh                # re-run the same command to resume after an interruption
```

## 3. Lawrencium: verify checksums, unpack and compare with the IGI inventories (one job)

```bash
TARGET=$HOME/IS110_ecology/Cycle_from_IGI
sbatch -A <account> -p <partition> -q <qos> -N 1 -c 12 --mem=64G -t 12:00:00 \
  -o ~/logs/verify_restore_%j.out \
  --wrap "bash $TARGET/transfer/migration_20261117/verify_on_target.sh $TARGET/transfer --restore $TARGET/data"
```

`data/` must not exist beforehand. The log must end with `ALL CHECKS PASSED` (about 3.5 h on IGI);
reports are written to `data/_verification/`.

## 4. Lawrencium: re-point the two links that still lead to IGI

```bash
cd $TARGET/data/IS110_downstream/cargo_is110
ln -sfn $TARGET/data/IS110_downstream__external/is110_analysis/intermediate intermediate
ln -sfn $TARGET/data/IS110_downstream__external/is110_analysis/src_data src_data
```

## 5. Lawrencium: code and environment

```bash
git clone https://github.com/KuangHu/nanopore_insertion_callout.git $TARGET/code/Cycle
cd $TARGET/code/Cycle && git checkout migration/lawrencium
mamba env create -n opfi -f envs/opfi.environment.yml      # another -n if an opfi env already exists
mamba install -n opfi -c bioconda samtools=1.17
cp config/lawrencium.env.example config/site.env
sinfo -p <partition> -o "%c %m %f" | sort -u                # cores / memory per node
```

In `config/site.env` set `CYCLE_ENV_PREFIX` to the environment's path (`conda env list` shows it),
`CYCLE_DATA_ROOT=$TARGET/data/IS_cycle`,
`CYCLE_HMM_DIR=$TARGET/data/assets/hmm`, `CYCLE_ECOR_DIR=$TARGET/data/ECOR_batch`, the Slurm
account, partition and QOS, and cores/memory from the `sinfo` line. Any required value left as
`CHANGE_ME` stops every job script with a message naming it.

## 6. Lawrencium: test

```bash
conda activate opfi
python -m pytest tests                                         # 42 passed on IGI
CYCLE_DRY_RUN=1 bin/cycle-sbatch scripts/is110_build_json.sh   # inspect the generated sbatch line
```

Then one single-sample pipeline run (`docs/lawrence_setup.md`, section 5). When steps 3-6 pass,
milestone 2 (restored and usable on Lawrencium) is reached; only then may the IGI sources be
considered for cleanup.
