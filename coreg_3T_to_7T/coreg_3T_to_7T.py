#!/usr/bin/env python3

from pathlib import Path
import os

# Define FSL version
FSL_VERSION = "6.0.7.11"
FREESURFER_VERSION = "7.4.1"

# Define subject and session
subject = "sub-003"
session_3T = "ses-04"
session_7T = "ses-03" # reference (7T): ses-04 for sub-001 and sub-002, ses-03 for sub-003

include_alignStep = True


############## Define file paths for 3T and 7T MRI
pdw_3T_path = Path(f"/data/pt_02262/data/TH_bids/testdata_Taechang/dcm_imported/{subject}/{session_3T}/anat")
pdw_3T_pattern = f"{subject}_{session_3T}_acq-PDw*echo*1_*part-mag*.nii"
pdw_3T_files = list(pdw_3T_path.glob(pdw_3T_pattern))
if not pdw_3T_files:
    raise FileNotFoundError(f"No file matching pattern {pdw_3T_pattern} found in {pdw_3T_path}")

pdw_7T_path = Path(f"/data/pt_02262/data/TH_bids/testdata_Taechang/LORAKS/{subject}/{session_7T}/anat")
pdw_7T_pattern = f"{subject}_{session_7T}_acq-PDw_rec-loraksRsos*echo-01*part-mag*.nii"
# pdw_7T_path = Path(f"/data/pt_02262/data/TH_bids/bids/{subject}/{session_7T}/anat")
# pdw_7T_pattern = f"{subject}_{session_7T}_acq-PDw*echo-1_*.nii"
pdw_7T_files = list(pdw_7T_path.glob(pdw_7T_pattern))
if not pdw_7T_files:
    raise FileNotFoundError(f"No file matching pattern {pdw_7T_pattern} found in {pdw_7T_path}")
pdw_7T_file = pdw_7T_files[0]  # Assuming only one 7T file is present
pdw_7T = str(pdw_7T_file)


############# Define output directory and create it if it doesn't exist
output_dir = Path(f"/data/pt_02262/data/TH_bids/testdata_Taechang/dcm_imported/derivatives/coreg_3T_to_7T/{subject}/{session_3T}/intermediate")
output_dir.mkdir(parents=True, exist_ok=True)


############# Preprocess 7T file with SynthStrip
pdw_7T_synthstrip = output_dir / f"{pdw_7T_file.stem}_synthstrip.nii.gz"
synthstrip_7T_command = f"FREESURFER --version {FREESURFER_VERSION} mri_synthstrip --i {pdw_7T_file} --o {pdw_7T_synthstrip}"
try:
    result = os.system(synthstrip_7T_command)
    if result != 0:
        raise OSError(f"SynthStrip failed for 7T file with exit code {result}")
    print(f"--> SynthStrip preprocessing complete for 7T file. Output saved to {pdw_7T_synthstrip}")
except OSError as e:
    print(f"--> Error during SynthStrip preprocessing for 7T file: {e}")


for idx, pdw_3T_file in enumerate(pdw_3T_files, start=1):

    ############# Preprocess each 3T file with SynthStrip 
    pdw_3T = str(pdw_3T_file)
    pdw_3T_synthstrip = output_dir / f"{pdw_3T_file.stem}_synthstrip.nii.gz"
    synthstrip_3T_command = f"FREESURFER --version {FREESURFER_VERSION} mri_synthstrip --i {pdw_3T} --o {pdw_3T_synthstrip}"
    try:
        result = os.system(synthstrip_3T_command)
        if result != 0:
            raise OSError(f"SynthStrip failed for 3T file {pdw_3T} with exit code {result}")
        print(f"--> SynthStrip preprocessing complete for 3T file {pdw_3T}. Output saved to {pdw_3T_synthstrip}")
    except OSError as e:
        print(f"--> Error during SynthStrip preprocessing for 3T file {pdw_3T}: {e}")

    print("------------------------------")
    ############# Align 3T to 7T using FLIRT (applyxfm)
    if include_alignStep:
        pdw_3T_aligned = output_dir / f"{pdw_3T_file.stem}_aligned.nii.gz"
        applyxfm_command = f"FSL --version {FSL_VERSION} flirt -in {pdw_3T_synthstrip} -ref {pdw_7T_synthstrip} -out {pdw_3T_aligned} -applyxfm -usesqform"
        try:
            result = os.system(applyxfm_command)
            if result != 0:
                raise OSError(f"FLIRT -applyxfm failed for 3T file {pdw_3T_synthstrip} with exit code {result}")
            print(f"--> FLIRT -applyxfm complete for 3T file {pdw_3T_synthstrip}. Output saved to {pdw_3T_aligned}")
        except OSError as e:
            print(f"--> Error during FLIRT -applyxfm for 3T file {pdw_3T_synthstrip}: {e}")


    print("------------------------------")
    # Coregistration of 3T to 7T
    pdw_3T_preprocessed = pdw_3T_aligned if include_alignStep else pdw_3T_synthstrip

    coreg_output = output_dir / f"coreg_pdw_3T{session_3T}_to_7T{session_7T}_run-{idx}.nii.gz"
    coreg_matrix = output_dir / f"coreg_pdw_3T{session_3T}_to_7T{session_7T}_run-{idx}.mat"
    print(f"Coregistration of {pdw_3T_preprocessed} to {pdw_7T_synthstrip}...")

    flirt_command = f"""
    FSL --version {FSL_VERSION} \
    flirt -in {pdw_3T_preprocessed} \
        -ref {pdw_7T_synthstrip} \
        -out {coreg_output} \
        -omat {coreg_matrix} \
        -dof 12 \
        -interp sinc \
        -cost mutualinfo
    """
    try:
        result = os.system(flirt_command)
        if result != 0:
            raise OSError(f"Command failed with exit code {result}")
        print(f"--> Coregistration complete for {pdw_3T_preprocessed}. Output saved to {coreg_output}")
    except OSError as e:
        print(f"--> Error during coregistration for {pdw_3T_preprocessed}: {e}")

    print("------------------------------")
    print("------------------------------")