#!/usr/bin/env python3

import os
import subprocess

##################### TODO #####################
# finds all subjects and sessions in the directory
# when batching the jobs, rather run each subject/session separately for parallelization




labels_file = "/data/u_kuegler_software/miniforge3/envs/qsmxt/lib/python3.8/site-packages/qsmxt/aseg_labels.csv"

# os.system(f"qsmxt /path/to/bids/dir --premade 'gre' \
#            --subjects sub-1 --sessions ses-1 \
#            --runs run-1 run-2 run-3 \
#            --labels_file {labels_file} --auto_yes")


# subprocess.run([
#     "qsmxt",
#     "/path/to/bids/dir",
#     "--premade", "gre",
#     "--subjects", "sub-1",
#     "--sessions", "ses-1",
#     "--runs", "run-1", "run-2", "run-3",
#     "--labels_file", labels_file,
#     "--auto_yes"
# ], check=True)