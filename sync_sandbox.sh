#/bin/bash
INPUT=$MYSANDBOX/hugoj/reproducing_jiarongs_plots/
SUBDIR="postprocess/"
OUTPUT=a0_jiarong_simu/

# add -n for dry run
rsync -av --include='*/' --include='*.py' --include='*.c' --include='*.txt' --exclude='*' $INPUT $OUTPUT
