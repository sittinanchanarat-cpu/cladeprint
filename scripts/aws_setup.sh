#!/usr/bin/env bash
# Prepare a fresh Ubuntu 24.04 (x86_64) EC2 machine for the census.
# Run from inside the cloned repository; safe to re-run (each step skips what is already there).
#
#   git clone https://<token>@github.com/sittinanchanarat-cpu/cladeprint.git && cd cladeprint
#   bash scripts/aws_setup.sh            # SMOKE=1 also runs the census script on 2 test runs
#
# Steps: system packages -> micromamba + the cladeprint environment (envs/genotype.yml) -> Rust ->
# cargo build --release -> species panel (pilot/fetch_refs.sh) -> tool check -> optional smoke test.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."  # the repository root (git may not be installed yet)
log() { echo -e "\n== $*"; }
SUDO=$([[ $EUID -eq 0 ]] || echo sudo)  # root in a container (testing): no sudo

log "system packages"
$SUDO apt-get update -qq
$SUDO env DEBIAN_FRONTEND=noninteractive apt-get install -y -qq build-essential pkg-config git curl unzip tmux htop bzip2 gawk > /dev/null
# gawk: Ubuntu's default awk is mawk, which lacks and() (used by the species screen); installing gawk
# makes it /usr/bin/awk

log "micromamba"
export MAMBA_ROOT_PREFIX=$HOME/micromamba
if ! command -v micromamba > /dev/null && [[ ! -x $HOME/.local/bin/micromamba ]]; then
  mkdir -p $HOME/.local/bin
  curl -sSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C $HOME/.local bin/micromamba
fi
export PATH=$HOME/.local/bin:$PATH
grep -q 'MAMBA_ROOT_PREFIX' $HOME/.bashrc || {
  echo 'export PATH=$HOME/.local/bin:$HOME/.cargo/bin:$PATH' >> $HOME/.bashrc
  echo "export MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX" >> $HOME/.bashrc
}

log "conda environment 'cladeprint'"
if ! micromamba env list | grep -q '^ *cladeprint '; then
  micromamba create -y -q -f envs/genotype.yml
fi

log "Rust and cladeprint"
if ! command -v cargo > /dev/null && [[ ! -x $HOME/.cargo/bin/cargo ]]; then
  curl -sSf https://sh.rustup.rs | sh -s -- -y -q --profile minimal
fi
export PATH=$HOME/.cargo/bin:$PATH
cargo build --release -q
./target/release/cladeprint --version

log "species panels (7 genomes of the Cryptococcus complex; 6 Sporothrix)"
if [[ ! -s pilot/ref/panel/complex.nuclear.txt || ! -s pilot/ref/panel/sporothrix.nuclear.txt ]]; then
  micromamba run -n cladeprint bash pilot/fetch_refs.sh
fi

log "tool check"
missing=0
for t in kmc kmc_tools fastq-dump minimap2 samtools curl gzip python; do
  if micromamba run -n cladeprint which $t > /dev/null 2>&1; then
    printf '  %-11s %s\n' $t "$(micromamba run -n cladeprint which $t)"
  else
    echo "  MISSING $t"; missing=1
  fi
done
awk --version 2>/dev/null | grep -q 'GNU Awk' || { echo "  awk is not GNU awk (need and()): sudo update-alternatives --set awk /usr/bin/gawk"; missing=1; }
(( missing == 0 )) || { echo "tools missing: fix before the census"; exit 1; }
echo "  cores $(nproc), memory $(free -g | awk '/Mem:/{print $2}') GB, free disk $(df -h --output=avail . | tail -1)"

if [[ ${SMOKE:-0} == 1 ]]; then
  log "smoke test: census script on 2 test runs"
  head -3 census/test_runs.tsv > census/smoke_runs.tsv
  OUT=smoke_out SHEET=smoke_runs.tsv JOBS=2 micromamba run -n cladeprint bash census/fetch_count.sh
  cat census/smoke_out/logs/timing.tsv
  rm -rf census/smoke_out census/smoke_runs.tsv
fi
log "ready. Next: tmux, then run the census scripts"
