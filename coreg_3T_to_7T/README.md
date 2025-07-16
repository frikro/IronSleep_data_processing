# Coregistration of the 3T TH data to the 7T data

Robust pipeline for coregistering 3T to 7T MRI data, including qMRI map transformation, using SLURM batch scripts.

- Coregistration of a 3T PDw echo-1 to the corresponding LORAKS-reconstructed echo-1 of a 7T session 
- Apply the calculated transformation to the 3T qMRI maps
- Reason: use MASSP segmentation of the 7T data for the 3T images

## Features:
- **Dual registration methods**: Support for both FSL FLIRT (affine only) and ANTS SyN (affine + nonlinear) registration
- **Auto-discovery**: Automatically discovers subjects and sessions from BIDS-like directory structures  
- **Robust error handling**: Comprehensive validation and error reporting
- **Batch processing**: SLURM-based batch job submission for multiple subjects/sessions
- **qMRI integration**: Automatic application of transforms to all qMRI maps after coregistration

## Usage:

### Quick Start (Auto-discovery):
```bash
# Auto-discover all subjects and sessions, use FLIRT (default)
./submit_coreg_batch.sh \
    /data/pt_02262/data/TH_bids/bids/derivatives/LCPCA_distCorr \
    /data/pt_02262/data/TH_bids/bids/derivatives/LORAKS_LCPCA_distCorr \
    /data/pt_02262/data/TH_bids/bids/derivatives/qMRI \
    /data/pt_02262/data/TH_bids/bids/derivatives/coreg_3T_to_7T

# Auto-discover all subjects and sessions, use ANTS SyN registration
./submit_coreg_batch.sh --use-ants \
    /data/pt_02262/data/TH_bids/bids/derivatives/LCPCA_distCorr \
    /data/pt_02262/data/TH_bids/bids/derivatives/LORAKS_LCPCA_distCorr \
    /data/pt_02262/data/TH_bids/bids/derivatives/qMRI \
    /data/pt_02262/data/TH_bids/bids/derivatives/coreg_3T_to_7T
```

### Specific Subjects/Sessions:
```bash
# Process specific subjects and sessions with ANTS
./submit_coreg_batch.sh --use-ants \
    -sub "sub-001,sub-002" \
    -ses3T "ses-05,ses-06" \
    -refSes "ses-04" \
    /data/pdw_3T /data/pdw_7T /data/qmri_3T /data/output

# Skip initial alignment step
./submit_coreg_batch.sh --no-align --use-ants \
    -sub "sub-001" -ses3T "ses-05" -refSes "ses-04" \
    /data/pdw_3T /data/pdw_7T /data/qmri_3T /data/output
```

### Options:
- `--use-ants`: Use ANTS SyN (affine + nonlinear) instead of FLIRT (affine only)
- `--no-align`: Skip initial FLIRT alignment step
- `--dry-run`: Show commands without submitting jobs
- `-sub`: Comma-separated list of subjects (auto-discovered if not specified)
- `-ses3T`: Comma-separated list of 3T sessions (auto-discovered if not specified)  
- `-refSes`: 7T reference session (defaults: sub-001/002→ses-04, sub-003→ses-03)
- `-t`: Delay between job submissions in seconds (default: 1)

## Registration Methods:

### FLIRT (Default - Affine Only):
- Fast, robust affine registration using mutual information
- Good for global alignment between 3T and 7T
- Suitable when anatomical differences are minimal

### ANTS SyN (Affine + Nonlinear):
- Multi-stage registration: Rigid → Affine → SyN nonlinear
- Better handles local anatomical differences
- More computationally intensive but higher accuracy
- Recommended for cases with significant B0 distortions or anatomical differences

## Output Structure:
```
output_directory/
├── sub-XXX/
│   └── ses-XX/
│       ├── intermediate/           # Coregistration outputs
│       │   ├── *_synthstrip.nii.gz      # Brain-extracted images
│       │   ├── coreg_pdw_*.nii.gz        # Coregistered PDw
│       │   ├── coreg_pdw_*.mat           # FLIRT transform matrix
│       │   ├── ants_*Warped.nii.gz       # ANTS warped image
│       │   ├── ants_*0GenericAffine.mat  # ANTS affine transform
│       │   ├── ants_*1Warp.nii.gz        # ANTS nonlinear warp
│       │   └── ants_*1InverseWarp.nii.gz # ANTS inverse warp
│       └── *_coreg.nii.gz          # Final qMRI maps in 7T space
```

## Pipeline Steps:
1. **Preprocessing**: Brain extraction using FreeSurfer SynthStrip
2. **Optional Alignment**: Initial FLIRT alignment using image headers (`-usesqform`)
3. **Registration**: 
   - FLIRT: 12-DOF affine with mutual information
   - ANTS: Rigid → Affine → SyN nonlinear with mutual information + cross-correlation
4. **qMRI Application**: Apply transforms to R1map, R2starmap, MTsat, PDmap

## Requirements:
- FSL 6.0.6 (via SCWRAP)
- FreeSurfer 7.4.1 (via SCWRAP) 
- ANTS 2.6.0 (via SCWRAP)
- SLURM environment

## Legacy Files (for reference):
- `coreg_3T_to_7T.py`: Original Python implementation
- `apply_coreg_qMRI.py`: Original Python qMRI application script




# Initial Tests:

```
./submit_coreg_batch.sh -sub sub-002 -ses3T ses-06,ses-07 -refSes ses-04 /data/pt_02262/data/TH_bids/bids/derivatives/LCPCA_distCorr /data/pt_02262/data/TH_bids/bids/derivatives/LORAKS_LCPCA_distCorr /data/pt_02262/data/TH_bids/bids/derivatives/qMRI /data/pt_02262/data/TH_bids/bids/derivatives/coreg_3T_to_7T/automated_flirt
```