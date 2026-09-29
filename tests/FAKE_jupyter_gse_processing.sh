#!/usr/bin/env bash
# FAKE stand-in for the lab's jupyter_gse_processing, used ONLY by the local tests.
# The real agent runs /booleanfs2/sahoo/Data/BooleanLab/Training/GSE_Processing_scripts/jupyter_gse_processing.
# Called the way the notebooks call it: bash <script> <folder name>, from inside that folder.
set -e
exec python3 "$(dirname "${BASH_SOURCE[0]}")/fake_support_files.py" "$1"
