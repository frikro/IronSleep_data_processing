# Coregistration of the 3T TH data to the 7T data

- coregistration of a 3T PDw echo-1 to the corresponding LORAKS-reconstructed echo-1 of a 7T session 
- apply the calculated warp field to the 3T qMRI maps
    - reason: use MASSP segmentation of the 7T data for the 3T images


## Usage:
- first run `coreg_3T_to_7T.py`
    - specify the correct subject and session
- after that, run `apply_coreg_qMRI.py` 
    - for the same session