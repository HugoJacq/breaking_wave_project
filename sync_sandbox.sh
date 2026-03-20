#/bin/bash
#
# add -n for dry run
#
# run with ./sync_sandbox.sh

# Jiarong sim
INPUT=$MYSANDBOX/hugoj/reproducing_jiarongs_plots/
SUBDIR="postprocess/"
OUTPUT=a0_jiarong_simu/
rsync -av --include='*/' --include='*.py' --include='*.c' --include='*.txt' --exclude='*' $INPUT $OUTPUT
rsync -av --include='*/' --include='*.py' --include='*.c' --include='*.txt' --exclude='*' $OUTPUT $INPUT

# verif_spectrum
INPUT="verif_spectrum/"
OUTPUT=$MYSANDBOX/hugoj/verif_spectrum/
rsync -av --include='*/' --include='*.py' --include='*.md' --include='*.txt' --include='Makefile' --exclude='*' $INPUT $OUTPUT
rsync -av --include='*/' --include='*.py' --include='*.md' --include='*.txt' --include='Makefile' --exclude='*' $OUTPUT $INPUT
