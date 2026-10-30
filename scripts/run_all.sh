#!/usr/bin/env bash
# Builds and runs every case (Re x mesh) with the official OpenFOAM v2406 Docker image.
# Usage: bash scripts/run_all.sh [parallel_jobs]   (run from the repository root)
set -euo pipefail
JOBS="${1:-8}"
IMAGE="opencfd/openfoam-run:2406"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Re  nu      endTime  writeInterval  meshes
CASES=(
  "100  0.01   30  10  32 64 128"
  "1000 0.001  80  10  32 64 128 256"
)

mkdir -p "$ROOT/runs"
list="$ROOT/runs/cases.txt"
: > "$list"
for spec in "${CASES[@]}"; do
  read -r re nu end wi meshes <<< "$spec"
  for n in $meshes; do
    dir="runs/Re${re}_N${n}"
    rm -rf "$ROOT/$dir"
    cp -r "$ROOT/case-template" "$ROOT/$dir"
    dt=$(python3 -c "print(0.5/$n)")   # Courant number 0.5 with lid speed 1 m/s
    sed -i.bak "s/__N__/$n/" "$ROOT/$dir/system/blockMeshDict"
    sed -i.bak -e "s/__END__/$end/" -e "s/__DT__/$dt/" -e "s/__WI__/$wi/" "$ROOT/$dir/system/controlDict"
    sed -i.bak "s/__NU__/$nu/" "$ROOT/$dir/constant/transportProperties"
    find "$ROOT/$dir" -name '*.bak' -delete
    echo "$dir" >> "$list"
  done
done

# The image entrypoint changes to the container home, so the case directory is entered explicitly.
# Meshes up to 128x128 run serially; 256x256 runs on 8 MPI processes (system/decomposeParDict).
grep -v '_N256$' "$list" | xargs -P "$JOBS" -I{} docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT:/work" "$IMAGE" \
  bash -c 'cd /work/{} && blockMesh > log.blockMesh 2>&1 && checkMesh > log.checkMesh 2>&1 && icoFoam > log.icoFoam 2>&1 && echo "done {}"'
grep '_N256$' "$list" | while read -r dir; do
  docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT:/work" "$IMAGE" bash -c "cd /work/$dir && blockMesh > log.blockMesh 2>&1 \
    && checkMesh > log.checkMesh 2>&1 && decomposePar > log.decomposePar 2>&1 \
    && mpirun -np 8 --oversubscribe icoFoam -parallel > log.icoFoam 2>&1 && echo done $dir"
done
