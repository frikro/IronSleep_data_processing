#!/bin/bash

# Script to cycle through subjects/sessions and submit SLURM jobs for 3T to 7T coregistration

usage() {
echo \
"
$(basename $0): Automatically submits SLURM jobs for 3T to 7T coregistration processing for specified subject/session combinations.

USAGE:
    $(basename $0) [options] <pdw_3T_directory> <pdw_7T_directory> <qmri_3T_directory> <output_directory>

OPTIONS:
    -h | --help: print help text and exit
    -sub SUBJECTS | --subjects SUBJECTS: comma-separated list (no spaces!) of subjects to process (e.g., sub-001,sub-002)
                                        If not specified, all subjects found in qmri_3T_directory will be processed
    -ses3T SESSIONS | --sessions-3T SESSIONS: comma-separated list (no spaces!) of 3T sessions to coregister to 7T space (e.g., ses-05,ses-06)
                                             If not specified, all 3T sessions found for each subject will be processed
    -refSes REFERENCE_SESSION | --reference-session REFERENCE_SESSION: 7T sessions to use as reference (e.g., ses-04)
                                                                        If not specified, default reference sessions will be used:
                                                                        sub-001, sub-002: ses-04; sub-003: ses-03, sub-005: ses-05
    -t SECONDS | --delay SECONDS: delay between job submissions in seconds (default: 1)
    --no-align: skip the initial alignment step (FLIRT applyxfm with usesqform)
    --use-ants: use ANTS SyN (affine + nonlinear) registration instead of FLIRT (affine only)
    --dry-run: show commands that would be executed without actually submitting jobs
    --qsm-dir DIRECTORY: (optional) Directory containing 3T QSM files to coregister
    --r2-dir DIRECTORY: (optional) Directory containing 3T R2 files to coregister
    --r2prime-dir DIRECTORY: (optional) Directory containing 3T R2' (R2prime) files to coregister

ARGUMENTS:
    pdw_3T_directory: Directory containing 3T PDw files (to calculate coregistration)
    pdw_7T_directory: Directory containing 7T PDw files (to calculate coregistration)
    qmri_3T_directory: Directory containing 3T qMRI files (to apply transforms)
    output_directory: Output directory for processed results

DESCRIPTION:
    The script submits SLURM jobs for coregistering 3T to 7T MRI data using FLIRT or ANTS.
    
    If no subjects or sessions are specified, the script will automatically discover all available
    subjects and 3T sessions by scanning the qmri_3T_directory for BIDS-structured data.
    
    For each subject/session combination, it:
    1. Preprocesses both 3T and 7T files with SynthStrip
    2. Optionally applies initial alignment using FLIRT -applyxfm -usesqform
    3. Performs final coregistration using either:
       - FLIRT with mutual information cost function (affine only, default)
       - ANTS SyN registration (affine + nonlinear, if --use-ants flag is used)
    4. Applies the calculated transformation to all qMRI files
    5. If specified, also applies the transformation to QSM, R2, and/or R2' files
    
    Creates output structure: output/sub-xxx/ses-xx/intermediate/ and output/sub-xxx/ses-xx/

    Expected input directory structures (BIDS-like):
      qMRI:    <qmri_dir>/<subject>/<session>/anat/<subject>_<session>_{R1map,R2starmap,MTsat,PDmap}.nii*
      QSM:     <qsm_dir>/<subject>/<session>/anat/coreg_toPDw/<subject>_<session>_mean_Chimap.nii*
      R2:      <r2_dir>/<subject>/<session>/anat/<subject>_<session>_*R2map.nii*
      R2prime: <r2prime_dir>/<subject>/<session>/anat/<subject>_<session>_*R2primemap.nii*

EXAMPLES:
    # Auto-discover all subjects and sessions
    $(basename $0) \\
        /data/pt_02262/data/TH_bids/bids/derivatives/LCPCA_distCorr \\
        /data/pt_02262/data/TH_bids/bids/derivatives/LORAKS_LCPCA_distCorr \\
        /data/pt_02262/data/TH_bids/bids/derivatives/qMRI \\
        /data/pt_02262/data/TH_bids/bids/derivatives/coreg_3T_to_7T

    # Specify particular subjects and sessions
    $(basename $0) -sub \"sub-001,sub-002\" -ses3T \"ses-05,ses-06\" -refSes \"ses-04\" \\
        /data/pdw_3T /data/pdw_7T /data/qmri_3T /data/output
    
    $(basename $0) -sub \"sub-001\" -ses3T \"ses-05\" -refSes \"ses-04\" \\
        /data/pt_02262/data/TH_bids/bids/derivatives/LCPCA_distCorr \\
        /data/pt_02262/data/TH_bids/bids/derivatives/LORAKS_LCPCA_distCorr \\
        /data/pt_02262/data/TH_bids/bids/derivatives/qMRI \\
        /data/pt_02262/data/TH_bids/bids/derivatives/coreg_3T_to_7T
    
    $(basename $0) --dry-run -sub \"sub-003\" -ses3T \"ses-05\" -refSes \"ses-03\" \\
        /data/pdw_3T /data/pdw_7T /data/qmri_3T /data/output
    
    $(basename $0) --no-align -sub \"sub-001\" -ses3T \"ses-05,ses-06\" -refSes \"ses-04\" \\
        /data/pdw_3T /data/pdw_7T /data/qmri_3T /data/output
    
    # Use ANTS SyN registration instead of FLIRT
    $(basename $0) --use-ants -sub \"sub-001\" -ses3T \"ses-05\" -refSes \"ses-04\" \\
        /data/pdw_3T /data/pdw_7T /data/qmri_3T /data/output
    
    # Also coregister QSM, R2 and R2' data
    $(basename $0) -sub \"sub-001\" -ses3T \"ses-05\" -refSes \"ses-04\" \\
        --qsm-dir /data/derivatives/QSM \\
        --r2-dir /data/derivatives/R2 \\
        --r2prime-dir /data/derivatives/R2prime \\
        /data/pdw_3T /data/pdw_7T /data/qmri_3T /data/output

AUTHOR:
    Niklas Kuegler (kuegler@cbs.mpg.de)
"
}

# Default parameters
delay=1
dry_run=false
include_align_step=true
use_ants=false
pdw_3T_directory=""
pdw_7T_directory=""
qmri_3T_directory=""
output_dir=""
subjects=""
sessions_3T=""
reference_session=""
qsm_directory=""
r2_directory=""
r2prime_directory=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        -sub|--subjects)
            subjects="$2"
            shift 2
            ;;
        -ses3T|--sessions-3T)
            sessions_3T="$2"
            shift 2
            ;;
        -refSes|--reference-session)
            reference_session="$2"
            shift 2
            ;;
        -t|--delay)
            delay="$2"
            shift 2
            ;;
        --no-align)
            include_align_step=false
            shift
            ;;
        --use-ants)
            use_ants=true
            shift
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        --qsm-dir)
            qsm_directory="$2"
            shift 2
            ;;
        --r2-dir)
            r2_directory="$2"
            shift 2
            ;;
        --r2prime-dir)
            r2prime_directory="$2"
            shift 2
            ;;
        -*)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
        *)
            if [[ -z "$pdw_3T_directory" ]]; then
                pdw_3T_directory="$1"
                shift
            elif [[ -z "$pdw_7T_directory" ]]; then
                pdw_7T_directory="$1"
                shift
            elif [[ -z "$qmri_3T_directory" ]]; then
                qmri_3T_directory="$1"
                shift
            elif [[ -z "$output_dir" ]]; then
                output_dir="$1"
                shift
            else
                echo "Error: Too many arguments specified"
                usage
                exit 1
            fi
            ;;
    esac
done

# Auto-adjust align step based on ANTS usage
if [[ "$use_ants" == "true" ]]; then
    include_align_step=false
    echo "ANTS registration selected - initial alignment step disabled"
fi

# Validation
if [[ -z "$pdw_3T_directory" ]]; then
    echo "Error: 3T PDw directory must be specified"
    usage
    exit 1
fi

if [[ ! -d "$pdw_3T_directory" ]]; then
    echo "Error: 3T PDw directory does not exist: $pdw_3T_directory"
    exit 1
fi

if [[ -z "$pdw_7T_directory" ]]; then
    echo "Error: 7T PDw directory must be specified"
    usage
    exit 1
fi

if [[ ! -d "$pdw_7T_directory" ]]; then
    echo "Error: 7T PDw directory does not exist: $pdw_7T_directory"
    exit 1
fi

if [[ -z "$qmri_3T_directory" ]]; then
    echo "Error: 3T qMRI directory must be specified"
    usage
    exit 1
fi

if [[ ! -d "$qmri_3T_directory" ]]; then
    echo "Error: 3T qMRI directory does not exist: $qmri_3T_directory"
    exit 1
fi

if [[ -z "$output_dir" ]]; then
    echo "Error: Output directory must be specified"
    usage
    exit 1
fi

if [[ -n "$qsm_directory" ]] && [[ ! -d "$qsm_directory" ]]; then
    echo "Error: QSM directory does not exist: $qsm_directory"
    exit 1
fi

if [[ -n "$r2_directory" ]] && [[ ! -d "$r2_directory" ]]; then
    echo "Error: R2 directory does not exist: $r2_directory"
    exit 1
fi

if [[ -n "$r2prime_directory" ]] && [[ ! -d "$r2prime_directory" ]]; then
    echo "Error: R2prime directory does not exist: $r2prime_directory"
    exit 1
fi

# Function to get default 7T reference session for a subject
get_default_reference_session() {
    local subject="$1"
    case "$subject" in
        "sub-001"|"sub-002")
            echo "ses-04"
            ;;
        "sub-003")
            echo "ses-03"
            ;;
        "sub-005")
            echo "ses-05" 
            ;;
        *)
            echo "Error: No default reference session defined for $subject. Please specify a reference session with -refSes option." >&2
            exit 1
            ;;
    esac
}

# Auto-discover subjects if not specified
if [[ -z "$subjects" ]]; then
    echo "No subjects specified. Auto-discovering subjects from: $qmri_3T_directory"
    discovered_subjects=()
    while IFS= read -r -d '' subject_dir; do
        subject=$(basename "$subject_dir")
        if [[ "$subject" =~ ^sub-[0-9]+$ ]]; then
            discovered_subjects+=("$subject")
        fi
    done < <(find "$qmri_3T_directory" -mindepth 1 -maxdepth 1 -type d -name "sub-*" -print0 2>/dev/null)
    
    if [[ ${#discovered_subjects[@]} -eq 0 ]]; then
        echo "Error: No subjects found in $qmri_3T_directory"
        exit 1
    fi
    
    # Sort subjects for consistent ordering
    IFS=$'\n' discovered_subjects=($(sort <<<"${discovered_subjects[*]}"))
    unset IFS
    
    subjects=$(IFS=','; echo "${discovered_subjects[*]}")
    echo "Found subjects: $subjects"
fi

# Convert comma-separated subjects to array
IFS=',' read -ra subject_array <<< "$subjects"

# Auto-discover sessions if not specified
if [[ -z "$sessions_3T" ]]; then
    echo "No 3T sessions specified. Auto-discovering sessions for each subject..."
    all_sessions=()
    
    for subject in "${subject_array[@]}"; do
        subject_sessions=()
        subject_qmri_dir="$qmri_3T_directory/$subject"
        
        if [[ -d "$subject_qmri_dir" ]]; then
            while IFS= read -r -d '' session_dir; do
                session=$(basename "$session_dir")
                if [[ "$session" =~ ^ses-[0-9]+$ ]]; then
                    # Check if this session has qMRI files
                    anat_dir="$session_dir/anat"
                    if [[ -d "$anat_dir" ]]; then
                        qmri_files=($(find "$anat_dir" -name "${subject}_${session}_*map.nii*" -type f 2>/dev/null))
                        if [[ ${#qmri_files[@]} -gt 0 ]]; then
                            subject_sessions+=("$session")
                        fi
                    fi
                fi
            done < <(find "$subject_qmri_dir" -mindepth 1 -maxdepth 1 -type d -name "ses-*" -print0 2>/dev/null)
            
            # Sort sessions for consistent ordering
            if [[ ${#subject_sessions[@]} -gt 0 ]]; then
                IFS=$'\n' subject_sessions=($(sort <<<"${subject_sessions[*]}"))
                unset IFS
                all_sessions+=("${subject_sessions[@]}")
                echo "  $subject: ${subject_sessions[*]}"
            else
                echo "  WARNING: No valid 3T sessions with qMRI files found for $subject"
            fi
        else
            echo "  WARNING: Subject directory not found: $subject_qmri_dir"
        fi
    done
    
    if [[ ${#all_sessions[@]} -eq 0 ]]; then
        echo "Error: No valid 3T sessions found for any subject"
        exit 1
    fi
    
    # Remove duplicates and sort
    IFS=$'\n' unique_sessions=($(printf '%s\n' "${all_sessions[@]}" | sort -u))
    unset IFS
    
    sessions_3T=$(IFS=','; echo "${unique_sessions[*]}")
    echo "Found 3T sessions: $sessions_3T"
fi

# Convert comma-separated sessions to array  
IFS=',' read -ra session_3T_array <<< "$sessions_3T"

# Create output directory if it doesn't exist
if [[ ! -d "$output_dir" ]]; then
    echo "Creating output directory: $output_dir"
    if [[ "$dry_run" == "false" ]]; then
        mkdir -p "$output_dir"
    fi
fi

# Create logs directory
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log_dir="$script_dir/logs"
if [[ "$dry_run" == "false" ]]; then
    mkdir -p "$log_dir"
fi

# Get the absolute path of the slurm script
slurm_script="$script_dir/coreg_3T_to_7T_slurm.sh"

if [[ ! -f "$slurm_script" ]]; then
    echo "Error: SLURM script not found at $slurm_script"
    exit 1
fi

echo "Subjects to process: ${subject_array[*]}"
echo "3T sessions to process: ${session_3T_array[*]}"
if [[ -n "$reference_session" ]]; then
    echo "7T reference session (specified): $reference_session"
else
    echo "7T reference session: Using defaults per subject (sub-001/sub-002: ses-04, sub-003: ses-03, sub-005: ses-05)"
fi
echo "Include alignment step: $include_align_step"
echo "Use ANTS registration: $use_ants"
echo "Output directory: $output_dir"
if [[ -n "$qsm_directory" ]]; then echo "QSM directory: $qsm_directory"; fi
if [[ -n "$r2_directory" ]]; then echo "R2 directory: $r2_directory"; fi
if [[ -n "$r2prime_directory" ]]; then echo "R2prime directory: $r2prime_directory"; fi

echo "=========================================="
echo "Jobs to be submitted:"

# Calculate and display job combinations
job_combinations=()
for subject in "${subject_array[@]}"; do
    for i in "${!session_3T_array[@]}"; do
        session_3T="${session_3T_array[$i]}"
        
        # Determine reference session
        if [[ -n "$reference_session" ]]; then
            ref_session="$reference_session"
        else
            ref_session=$(get_default_reference_session "$subject")
        fi
        
        # Only add combination if subject has this session
        subject_qmri_dir="$qmri_3T_directory/$subject/$session_3T/anat"
        if [[ -d "$subject_qmri_dir" ]]; then
            qmri_files=($(find "$subject_qmri_dir" -name "${subject}_${session_3T}_*map.nii*" -type f 2>/dev/null))
            if [[ ${#qmri_files[@]} -gt 0 ]]; then
                job_combinations+=("$subject:$session_3T:$ref_session")
                echo "$subject $session_3T -> $ref_session"
            fi
        fi
    done
done
echo "=========================================="

if [[ ${#job_combinations[@]} -eq 0 ]]; then
    echo "Error: No valid job combinations found"
    exit 1
fi

# Submit jobs
job_counter=1
total_jobs=${#job_combinations[@]}

for job_combo in "${job_combinations[@]}"; do
    IFS=':' read -ra combo_parts <<< "$job_combo"
    subject="${combo_parts[0]}"
    session_3T="${combo_parts[1]}"
    ref_session="${combo_parts[2]}"
    
    echo ""
    echo "Submitting job $job_counter/$total_jobs: $subject $session_3T -> $ref_session"
    
    # Check if output already exists
    expected_output_pattern="${output_dir}/${subject}/${session_3T}/intermediate/coreg_pdw_3T${session_3T}_to_7T${ref_session}.nii.gz"
    existing_outputs=($(ls $expected_output_pattern 2>/dev/null))
    
    if [[ ${#existing_outputs[@]} -gt 0 ]]; then
        echo "  INFO: Output files already exist for $subject $session_3T -> $ref_session. Skipping."
        ((job_counter++))
        continue
    fi
    
    # Prepare SLURM command
    SLURM_PARTITIONS="short,group_servers,gr_weiskopf"
    slurm_cmd="sbatch -p \"$SLURM_PARTITIONS\" \"$slurm_script\" \"$subject\" \"$session_3T\" \"$ref_session\" \"$pdw_3T_directory\" \"$pdw_7T_directory\" \"$qmri_3T_directory\" \"$output_dir\" \"$include_align_step\" \"$use_ants\" \"$qsm_directory\" \"$r2_directory\" \"$r2prime_directory\""
    
    if [[ "$dry_run" == "false" ]]; then
        # Submit the job
        out=$(eval $slurm_cmd)
        echo "  $out"
        
        # Add delay between submissions (except for the last job)
        if [[ $job_counter -lt $total_jobs ]]; then
            sleep "$delay"
        fi
    else
        echo "  DRY RUN: Would submit job with command: $slurm_cmd"
    fi
    
    ((job_counter++))
done

echo ""
echo "=========================================="
echo "Batch submission completed!"
echo "Total jobs: $total_jobs"
if [[ "$dry_run" == "false" ]]; then
    echo "Check job status with: squeue -u \$USER"
    echo "Monitor logs in: $log_dir/"
else
    echo "This was a dry run - no jobs were actually submitted"
fi
echo "=========================================="
