#!/usr/bin/env python3

import os
import re


################## define variables ##################

source_dir_parent = "/data/pt_02262/data/TH_bids/testdata_Taechang/LORAKS"
temp_dir_parent = "/data/pt_02262/data/TH_bids/testdata_Taechang/LORAKS/derivatives/qsmxt"
subjects_list = ["001", "002", "004"] ### TODO: find all "sub-*" and find all "ses-*" in the directory
sessions_list = ["04"]

#### TODO ####
# typical naming: sub-1_acq-mygrea_run-1_echo-1_part-mag_MEGRE.nii
# rename the LORAKS reconstructed files to the typical naming (symlinks)
# place the temporary directory in the bids folder -> derivatives needs to be at the correct spot


################## functions ##################

def create_renamed_symlinks(session_dir, temp_dir, contrasts, suffix):
    """
    Create renamed symbolic links for files in a session directory.
    This function scans the specified session directory for files that match the given suffix and 
    contain "rec-loraks_" in their names. It then creates symbolic links to these files in a 
    temporary directory, renaming them to be compliant with QSMxT.
    
    Parameters:
        session_dir (str): The directory containing the original files.
        temp_dir (str): The directory where the symbolic links will be created.
        contrasts (str or list of str): The contrast(s) to look for in the filenames. If a single 
                                        contrast is provided as a string, it will be converted to a list.
        suffix (str): The suffix to look for in the filenames. If the suffix is not "MEGRE", it will 
                    be replaced by "MEGRE" in the new filenames.
    Returns:
        None
    """
    

    os.makedirs(temp_dir, exist_ok=True)

    for fname in os.listdir(session_dir):

        lower_fname = fname.casefold()
        if suffix.casefold() in lower_fname and "rec-loraks_" in lower_fname:
            
            if isinstance(contrasts, str):
                contrasts = [contrasts]
            
            for idx, contr in enumerate(contrasts):
                if contr.casefold() in lower_fname:
                    
                    ### necessary in older versions of qsmxt
                    # if suffix != "MEGRE":
                    #     # Replace suffix by "MEGRE"
                    #     new_fname = re.sub(suffix, "MEGRE", fname, flags=re.IGNORECASE)
                    # else:
                    #     new_fname = fname
                    new_fname = fname

                    # Insert "_run-*" after whichever contrast is found 
                    # (only if multiple contrasts are specified)
                    if len(contrasts) > 1:
                        pos = new_fname.casefold().find(contr.casefold())
                        if pos != -1:
                            new_fname = (
                                new_fname[:pos+len(contr)]
                                + f"_run-{idx+1}"
                                + new_fname[pos+len(contr):]
                            )

                    os.system(f"ln -s {os.path.join(session_dir, fname)} {os.path.join(temp_dir, new_fname)}")



################## main ##################

for subject in subjects_list:
    for session in sessions_list:

        session_dir = os.path.join(source_dir_parent, f"sub-{subject}", f"ses-{session}", "anat")

        temp_dir_mpm = os.path.join(temp_dir_parent, "temp_qsmxt", f"sub-{subject}", f"ses-{session}", "anat")
        suffix = "MPM"
        contrasts = ["T1w", "PDw", "MTw"] # renaming into run-1, run-2, run-3
        create_renamed_symlinks(session_dir, temp_dir_mpm, contrasts, suffix)


        temp_dir_slab = os.path.join(temp_dir_parent, "temp_qsmxt_slab", f"sub-{subject}", f"ses-{session}", "anat")
        suffix = "MEGRE"
        contrasts = "ernst" # no run if only one contrast
        create_renamed_symlinks(session_dir, temp_dir_slab, contrasts, suffix)

