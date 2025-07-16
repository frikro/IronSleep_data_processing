#!/bin/bash

#
#SBATCH -c 16					# 16 cores
#SBATCH --mem 64G				# estimated 32G RAM
#SBATCH --time 180				# estimated 90 minutes maximum
#SBATCH -o /data/u_kuegler_software/git/ironsleep_data_processing/coreg_3T_to_7T/logs/%j.out	# redirect the output
#

# Script for coregistering 3T to 7T MRI data using FLIRT or ANTS
# Arguments:
# $1: subject (e.g., sub-001)
# $2: session_3T (e.g., ses-05)
# $3: reference_session (e.g., ses-04)
# $4: 3Tpdw_directory (bids directory containing 3T Pdw files)
# $5: 7Tpdw_directory (bids directory containing 7T Pdw files)
# $6: 3Tqmri_directory (bids directory containing 3T qMRI files)
# $7: output_dir (base output directory)
# $8: include_align_step (true/false)
# $9: use_ants (true/false) - if true, use ANTS SyN instead of FLIRT

subject=$1
session_3T=$2
reference_session=$3
pdw_3T_directory=$4
pdw_7T_directory=$5
qmri_3T_directory=$6
output_dir=$7
include_align_step=$8
use_ants=$9

# Define FSL and FreeSurfer versions
FREESURFER_ENV="SCWRAP freesurfer 7.4.1"
FSL_ENV="SCWRAP fsl 6.0.6"
ANTS_ENV="SCWRAP ants 2.6.0"

echo "=== Coregistration Job Started ==="
echo "Subject: $subject"
echo "3T Session: $session_3T"
echo "7T reference Session: $reference_session"
echo "3T Pdw Directory: $pdw_3T_directory"
echo "7T Pdw Directory: $pdw_7T_directory"
echo "3T QMRI Directory: $qmri_3T_directory"
echo "Output directory: $output_dir"
echo "Include align step: $include_align_step"
echo "Use ANTS SyN: $use_ants"
echo "=================================="

# Define file paths for 3T and 7T MRI
# 3T Pdw files
pdw_3T_path="$pdw_3T_directory/${subject}/${session_3T}/anat"
pdw_3T_pattern="${subject}_${session_3T}_acq-PDw*echo-01_*part-mag*.nii"

# 7T Pdw files
pdw_7T_path="$pdw_7T_directory/${subject}/${reference_session}/anat"
pdw_7T_pattern="${subject}_${reference_session}_acq-PDw_rec-loraksRsos*echo-01*part-mag*.nii"

# Create output directory structure
intermediate_dir="${output_dir}/${subject}/${session_3T}/intermediate"
mkdir -p "$intermediate_dir"

echo "Searching for 3T files in: $pdw_3T_path"
echo "3T pattern: $pdw_3T_pattern"

# Find 3T files
pdw_3T_files=($(find "$pdw_3T_path" -name "$pdw_3T_pattern" 2>/dev/null))
if [[ ${#pdw_3T_files[@]} -eq 0 ]]; then
    echo "ERROR: No 3T files matching pattern $pdw_3T_pattern found in $pdw_3T_path"
    exit 1
fi

if [[ ${#pdw_3T_files[@]} -gt 1 ]]; then
    echo "ERROR: Found ${#pdw_3T_files[@]} 3T PDw files. This indicates multiple runs exist."
    echo "Found files:"
    for file in "${pdw_3T_files[@]}"; do
        echo "  $file"
    done
    echo ""
    echo "Multiple 3T PDw files detected (likely different run- entities)."
    echo "Cannot determine which specific run was used to create the qMRI maps."
    echo "Please specify which run to use by modifying the search pattern or"
    echo "ensure only one 3T PDw file exists for this session."
    echo ""
    echo "You may need to:"
    echo "1. Check which run was used during qMRI processing"
    echo "2. Update the script to specify the exact run (e.g., run-01)"
    echo "3. Or process each run separately with corresponding qMRI maps"
    exit 1
fi

pdw_3T_file="${pdw_3T_files[0]}"
echo "Found 3T file: $pdw_3T_file"

echo "Searching for 7T files in: $pdw_7T_path"
echo "7T pattern: $pdw_7T_pattern"

# Find 7T files
pdw_7T_files=($(find "$pdw_7T_path" -name "$pdw_7T_pattern" 2>/dev/null))
if [[ ${#pdw_7T_files[@]} -eq 0 ]]; then
    echo "ERROR: No 7T files matching pattern $pdw_7T_pattern found in $pdw_7T_path"
    exit 1
fi

if [[ ${#pdw_7T_files[@]} -gt 1 ]]; then
    echo "ERROR: Multiple 7T files found. Expected exactly one file."
    echo "Found files:"
    for file in "${pdw_7T_files[@]}"; do
        echo "  $file"
    done
    exit 1
fi

pdw_7T_file="${pdw_7T_files[0]}"
echo "Using 7T reference file: $pdw_7T_file"

# Extract basename for 7T file
pdw_7T_basename=$(basename "$pdw_7T_file" .nii.gz)
pdw_7T_basename=$(basename "$pdw_7T_basename" .nii)

# Preprocess 7T file with SynthStrip
echo "=== Preprocessing 7T file with SynthStrip ==="
pdw_7T_synthstrip="${intermediate_dir}/${pdw_7T_basename}_synthstrip.nii.gz"

if [[ ! -f "$pdw_7T_synthstrip" ]]; then
    echo "Running SynthStrip on 7T file..."
    $FREESURFER_ENV mri_synthstrip --i "$pdw_7T_file" --o "$pdw_7T_synthstrip"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: SynthStrip failed for 7T file"
        exit 1
    fi
    echo "SynthStrip preprocessing complete for 7T file. Output saved to $pdw_7T_synthstrip"
else
    echo "7T SynthStrip output already exists: $pdw_7T_synthstrip"
fi

# Process 3T file
echo ""
echo "=== Preprocessing 3T file with Synthstrip ==="

# Extract basename for 3T file
pdw_3T_basename=$(basename "$pdw_3T_file" .nii.gz)
pdw_3T_basename=$(basename "$pdw_3T_basename" .nii)

# Preprocess 3T file with SynthStrip
echo "Running SynthStrip on 3T file..."
pdw_3T_synthstrip="${intermediate_dir}/${pdw_3T_basename}_synthstrip.nii.gz"
    
if [[ ! -f "$pdw_3T_synthstrip" ]]; then
    $FREESURFER_ENV mri_synthstrip --i "$pdw_3T_file" --o "$pdw_3T_synthstrip"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: SynthStrip failed for 3T file $pdw_3T_file"
        continue
    fi
    echo "SynthStrip preprocessing complete for 3T file. Output saved to $pdw_3T_synthstrip"
else
    echo "3T SynthStrip output already exists: $pdw_3T_synthstrip"
fi
    
# Align 3T to 7T using FLIRT (applyxfm) if requested
echo "-------------------------------"
if [[ "$include_align_step" == "true" ]]; then
    echo "Running FLIRT alignment step..."
    pdw_3T_aligned="${intermediate_dir}/${pdw_3T_basename}_aligned.nii.gz"
    
    if [[ ! -f "$pdw_3T_aligned" ]]; then
        $FSL_ENV flirt -in "$pdw_3T_synthstrip" -ref "$pdw_7T_synthstrip" -out "$pdw_3T_aligned" -applyxfm -usesqform
        if [[ $? -ne 0 ]]; then
            echo "ERROR: FLIRT -applyxfm failed for 3T file $pdw_3T_synthstrip"
            continue
        fi
        echo "FLIRT -applyxfm complete. Output saved to $pdw_3T_aligned"
    else
        echo "3T aligned output already exists: $pdw_3T_aligned"
    fi
    pdw_3T_preprocessed="$pdw_3T_aligned"
else
    echo "Skipping alignment step"
    pdw_3T_preprocessed="$pdw_3T_synthstrip"
fi

echo "------------------------------"    
# Coregistration of 3T to 7T
if [[ "$use_ants" == "true" ]]; then
    echo "Running ANTS SyN coregistration (affine + nonlinear)..."
    
    # Define ANTS output files
    ants_prefix="${intermediate_dir}/ants_3T${session_3T}_to_7T${reference_session}_"
    ants_output="${ants_prefix}Warped.nii.gz"
    ants_affine="${ants_prefix}0GenericAffine.mat"
    ants_warp="${ants_prefix}1Warp.nii.gz"
    ants_inverse_warp="${ants_prefix}1InverseWarp.nii.gz"
    
    # For compatibility, create symlink with expected name
    coreg_output="${intermediate_dir}/coreg_pdw_3T${session_3T}_to_7T${reference_session}.nii.gz"
    
    if [[ ! -f "$ants_output" ]]; then
        echo "Running ANTS registration with SyN transformation..."
        echo "Fixed image (7T): $pdw_7T_synthstrip"
        echo "Moving image (3T): $pdw_3T_preprocessed"
        
        $ANTS_ENV antsRegistration \
            --dimensionality 3 \
            --float 0 \
            --output [$ants_prefix,$ants_output] \
            --interpolation Linear \
            --winsorize-image-intensities [0.005,0.995] \
            --use-histogram-matching 0 \
            --initial-moving-transform [$pdw_7T_synthstrip,$pdw_3T_preprocessed,1] \
            --transform Rigid[0.1] \
            --metric MI[$pdw_7T_synthstrip,$pdw_3T_preprocessed,1,32,Regular,0.25] \
            --convergence [1000x500x250x100,1e-6,10] \
            --shrink-factors 8x4x2x1 \
            --smoothing-sigmas 3x2x1x0vox \
            --transform Affine[0.1] \
            --metric MI[$pdw_7T_synthstrip,$pdw_3T_preprocessed,1,32,Regular,0.25] \
            --convergence [1000x500x250x100,1e-6,10] \
            --shrink-factors 8x4x2x1 \
            --smoothing-sigmas 3x2x1x0vox \
            --transform SyN[0.1,3,0] \
            --metric CC[$pdw_7T_synthstrip,$pdw_3T_preprocessed,1,4] \
            --convergence [100x70x50x20,1e-6,10] \
            --shrink-factors 8x4x2x1 \
            --smoothing-sigmas 3x2x1x0vox
        
        # other possible interpolation: Linear, BSPline, CosineWindowedSinc

        if [[ $? -ne 0 ]]; then
            echo "ERROR: ANTS registration failed for $pdw_3T_preprocessed"
            exit 1
        fi
        
        # Create symlink for compatibility
        ln -sf "$(basename "$ants_output")" "$coreg_output"
        
        echo "ANTS registration complete. Output saved to $ants_output"
        echo "Affine transform saved to $ants_affine"
        echo "Nonlinear warp saved to $ants_warp"
        echo "Inverse warp saved to $ants_inverse_warp"
    else
        echo "ANTS registration output already exists: $ants_output"
        # Ensure symlink exists
        if [[ ! -L "$coreg_output" ]]; then
            ln -sf "$(basename "$ants_output")" "$coreg_output"
        fi
    fi
    
    # Set variables for qMRI processing
    transform_type="ants"
    coreg_matrix=""  # ANTS doesn't use simple matrix
    
else
    echo "Running FLIRT coregistration (affine only)..."
    coreg_output="${intermediate_dir}/coreg_pdw_3T${session_3T}_to_7T${reference_session}.nii.gz"
    coreg_matrix="${intermediate_dir}/coreg_pdw_3T${session_3T}_to_7T${reference_session}.mat"
        
    if [[ ! -f "$coreg_output" ]]; then
        echo "Coregistering $pdw_3T_preprocessed to $pdw_7T_synthstrip..."
        
        $FSL_ENV flirt \
            -in "$pdw_3T_preprocessed" \
            -ref "$pdw_7T_synthstrip" \
            -out "$coreg_output" \
            -omat "$coreg_matrix" \
            -dof 12 \
            -interp sinc \
            -cost mutualinfo
            
        if [[ $? -ne 0 ]]; then
            echo "ERROR: FLIRT coregistration failed for $pdw_3T_preprocessed"
            exit 1
        fi
        echo "FLIRT coregistration complete. Output saved to $coreg_output"
        echo "Transformation matrix saved to $coreg_matrix"
    else
        echo "FLIRT coregistration output already exists: $coreg_output"
    fi
    
    # Set variables for qMRI processing
    transform_type="flirt"
fi
    
echo "------------------------------"

echo ""
echo "=== Coregistration Job Completed ==="
echo "All outputs saved in: $intermediate_dir"
echo "====================================="


echo ""
echo "=== Applying Co-registration to 3T qMRI files ==="

# Define input qMRI directory and file patterns
qmri_input_dir="${qmri_3T_directory}/${subject}/${session_3T}/anat"
qmri_output_dir="${output_dir}/${subject}/${session_3T}"
mkdir -p "$qmri_output_dir"

echo "Searching for qMRI files in: $qmri_input_dir"

# Define qMRI file patterns to process
qmri_patterns=(
    "${subject}_${session_3T}_R1map.nii*"
    "${subject}_${session_3T}_R2starmap.nii*"
    "${subject}_${session_3T}_MTsat.nii*"
    "${subject}_${session_3T}_PDmap.nii*"
)

# Find all qMRI files
qmri_files=()
for pattern in "${qmri_patterns[@]}"; do
    while IFS= read -r -d '' file; do
        qmri_files+=("$file")
    done < <(find "$qmri_input_dir" -name "$pattern" -print0 2>/dev/null)
done

if [[ ${#qmri_files[@]} -eq 0 ]]; then
    echo "WARNING: No qMRI files found matching patterns in $qmri_input_dir"
    echo "Patterns searched:"
    for pattern in "${qmri_patterns[@]}"; do
        echo "  $pattern"
    done
    echo "Skipping qMRI coregistration application."
else
    echo "Found ${#qmri_files[@]} qMRI files to process"
    
    # Check if coregistration transforms exist
    if [[ "$transform_type" == "ants" ]]; then
        if [[ ! -f "$ants_affine" ]] || [[ ! -f "$ants_warp" ]]; then
            echo "ERROR: ANTS transforms not found:"
            echo "  Affine: $ants_affine"
            echo "  Warp: $ants_warp"
            echo "Cannot apply coregistration to qMRI files."
            exit 1
        fi
        echo "Using ANTS transforms:"
        echo "  Affine: $ants_affine"
        echo "  Warp: $ants_warp"
    else
        if [[ ! -f "$coreg_matrix" ]]; then
            echo "ERROR: FLIRT coregistration matrix not found: $coreg_matrix"
            echo "Cannot apply coregistration to qMRI files."
            exit 1
        fi
        echo "Using FLIRT coregistration matrix: $coreg_matrix"
    fi
    
    echo "Using 7T reference file: $pdw_7T_synthstrip"
    
    # Process each qMRI file
    for qmri_file in "${qmri_files[@]}"; do
        echo ""
        echo "Processing qMRI file: $(basename "$qmri_file")"
        
        # Extract basename
        qmri_basename=$(basename "$qmri_file" .nii.gz)
        qmri_basename=$(basename "$qmri_basename" .nii)
        
        # Define output files
        qmri_aligned="${qmri_output_dir}/${qmri_basename}_aligned.nii.gz"
        qmri_coreg_output="${qmri_output_dir}/${qmri_basename}_coreg.nii.gz"
        
        # Apply initial alignment if requested
        if [[ "$include_align_step" == "true" ]]; then
            echo "  Applying initial alignment..."
            if [[ ! -f "$qmri_aligned" ]]; then
                $FSL_ENV flirt \
                    -in "$qmri_file" \
                    -ref "$pdw_7T_synthstrip" \
                    -out "$qmri_aligned" \
                    -applyxfm \
                    -usesqform
                
                if [[ $? -ne 0 ]]; then
                    echo "  ERROR: FLIRT -applyxfm failed for qMRI file $qmri_file"
                    continue
                fi
                echo "  Initial alignment complete. Output saved to $qmri_aligned"
            else
                echo "  Initial alignment output already exists: $qmri_aligned"
            fi
            qmri_preprocessed="$qmri_aligned"
        else
            echo "  Skipping initial alignment step"
            qmri_preprocessed="$qmri_file"
        fi
        
        # Apply coregistration transform
        echo "  Applying coregistration transform..."
        if [[ ! -f "$qmri_coreg_output" ]]; then
            if [[ "$transform_type" == "ants" ]]; then
                echo "  Using ANTS transform..."
                $ANTS_ENV antsApplyTransforms \
                    -d 3 \
                    -i "$qmri_preprocessed" \
                    -r "$pdw_7T_synthstrip" \
                    -o "$qmri_coreg_output" \
                    -t "$ants_warp" \
                    -t "$ants_affine" \
                    --interpolation Linear
                
                if [[ $? -ne 0 ]]; then
                    echo "  ERROR: ANTS transform application failed for $qmri_file"
                    continue
                fi
                echo "  ANTS transform applied successfully. Output saved to $qmri_coreg_output"
            else
                echo "  Using FLIRT transform..."
                $FSL_ENV flirt \
                    -in "$qmri_preprocessed" \
                    -ref "$pdw_7T_synthstrip" \
                    -applyxfm \
                    -init "$coreg_matrix" \
                    -out "$qmri_coreg_output"
                
                if [[ $? -ne 0 ]]; then
                    echo "  ERROR: FLIRT coregistration application failed for $qmri_file"
                    continue
                fi
                echo "  FLIRT transform applied successfully. Output saved to $qmri_coreg_output"
            fi
        else
            echo "  Coregistration output already exists: $qmri_coreg_output"
        fi
        
        echo "  ------------------------------"
    done
fi

echo ""
echo "=== qMRI Coregistration Application Completed ==="
echo "qMRI outputs saved in: $qmri_output_dir"
echo "=============================================="

echo ""
echo "=== FULL JOB COMPLETED ==="
echo "PDw coregistration outputs: $intermediate_dir"
echo "qMRI coregistration outputs: $qmri_output_dir"
echo "=========================="

