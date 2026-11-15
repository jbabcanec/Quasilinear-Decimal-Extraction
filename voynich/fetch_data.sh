#!/usr/bin/env bash
# Fetches the third-party data this study reads (not redistributed here).
#   - ZL3b transliteration (Zandbergen & Landini), word boxes (voynichese.com),
#     control corpora and the Naibbe / self-citation control texts,
#   as archived in Rozanova & Temerev's reproducibility repo, pinned to a commit.
set -euo pipefail
cd "$(dirname "$0")"
PIN=956a7c4fc39981f4d116fa3f4edfccce6d065571
if [ ! -d external/voynich-units/.git ]; then
  mkdir -p external
  GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/lrozanova/voynich-units external/voynich-units
fi
git -C external/voynich-units checkout -q "$PIN"
echo "data ready at external/voynich-units (pinned $PIN)"
