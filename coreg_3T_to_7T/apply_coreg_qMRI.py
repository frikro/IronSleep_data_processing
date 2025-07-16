#!/usr/bin/env python3

from pathlib import Path
import os

# Define FSL version
FSL_VERSION = "6.0.7.11"

# Define subject, session, and input files
subject = "sub-001"
session_3T = "ses-08"
session_7T = "ses-04" # reference (7T): ses-04 for sub-001 and sub-002, ses-03 for sub-003

include_alignStep = True

# Define input directory and files to process
# input_dir = Path(f"/data/pt_02262/data/TH_bids/testdata_Taechang/dcm_imported/derivatives/qMRI/{subject}/{session_3T}/anat")
input_dir = Path(f"/data/pt_02262/data/TH_bids/bids/derivatives/qMRI/{subject}/{session_3T}/anat")
input_files_patterns = [
    f"{subject}_{session_3T}_R1map.nii",
    f"{subject}_{session_3T}_R2starmap.nii",
    f"{subject}_{session_3T}_MTsat.nii",
    f"{subject}_{session_3T}_PDmap.nii",
    # "fit_r2.nii"
]

input_files = []
for pattern in input_files_patterns:
    input_files.extend(input_dir.glob(pattern))

if not input_files:
    raise FileNotFoundError(f"No files matching patterns {input_files_patterns} found in {input_dir}")

# Define coregistration matrix and output directory
# coreg_matrix_dir = Path(f"/data/pt_02262/data/TH_bids/testdata_Taechang/dcm_imported/derivatives/coreg_3T_to_7T/{subject}/{session_3T}/intermediate")
coreg_matrix_dir = Path(f"/data/pt_02262/data/TH_bids/bids/derivatives/coreg_3T_to_7T/flirt_affine/{subject}/{session_3T}/intermediate")
coreg_matrix_pattern = f"coreg_pdw_3T{session_3T}_to_7T{session_7T}_run-*.mat"
coreg_matrices = list(coreg_matrix_dir.glob(coreg_matrix_pattern))
if not coreg_matrices:
    raise FileNotFoundError(f"No coregistration matrices matching pattern {coreg_matrix_pattern} found in {coreg_matrix_dir}")

# Find reference file for coregistration
reference_file_7T_pattern = f"{subject}_{session_7T}*_acq-PDw_rec-loraksRsos_echo-01*_synthstrip.nii.gz"
reference_file_7T = list(coreg_matrix_dir.glob(reference_file_7T_pattern))
if not reference_file_7T:
    raise FileNotFoundError(f"No reference file matching pattern {reference_file_7T_pattern} found in {coreg_matrix_dir}")
reference_file_7T = str(reference_file_7T[0])  # Assuming only one reference file is present


# output_dir = Path(f"/data/pt_02262/data/TH_bids/testdata_Taechang/dcm_imported/derivatives/coreg_3T_to_7T/{subject}/{session_3T}")
output_dir = Path(f"/data/pt_02262/data/TH_bids/bids/derivatives/coreg_3T_to_7T/flirt_affine/{subject}/{session_3T}")
output_dir.mkdir(parents=True, exist_ok=True)

# Apply coregistration matrix to each input file
for idx, input_file in enumerate(input_files, start=1):
    coreg_matrix = coreg_matrices[-1]  # Use the last coregistration matrix
    coreg_output = output_dir / f"{input_file.stem}_coreg.nii.gz"


    ############# Align 3T to 7T using FLIRT (applyxfm)
    if include_alignStep:
        input_file_aligned = output_dir / f"{input_file.stem}_aligned.nii.gz"
        applyxfm_command = f"FSL --version {FSL_VERSION} flirt -in {input_file} -ref {reference_file_7T} -out {input_file_aligned} -applyxfm -usesqform"
        try:
            result = os.system(applyxfm_command)
            if result != 0:
                raise OSError(f"FLIRT -applyxfm failed for 3T file {input_file} with exit code {result}")
            print(f"--> FLIRT -applyxfm complete for 3T file {input_file}. Output saved to {input_file_aligned}")
        except OSError as e:
            print(f"--> Error during FLIRT -applyxfm for 3T file {input_file}: {e}")

    print("-----")
    print(f"Applying coregistration matrix {coreg_matrix} to {input_file}...")

    input_file_preprocessed = input_file_aligned if include_alignStep else input_file

    flirt_command = f"""
    FSL --version {FSL_VERSION} \
    flirt -in {input_file_preprocessed} \
        -ref {reference_file_7T} \
        -applyxfm -init {coreg_matrix} \
        -out {coreg_output}
    """

    try:
        result = os.system(flirt_command)
        if result != 0:
            raise OSError(f"FLIRT command failed for {input_file} with exit code {result}")
        print(f"--> Coregistration applied successfully. Output saved to {coreg_output}")
    except OSError as e:
        print(f"--> Error applying coregistration for {input_file}: {e}")

    print("------------------------------")
