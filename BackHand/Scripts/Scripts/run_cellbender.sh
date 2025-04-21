#!/bin/bash

# Arguments passed from Python script
donor=$1
base_data_dir=$2
expectedCell=$3  # Additional argument for expectedCell
pipeline_used=$4

if [ "$pipeline_used" == "count" ]; then
    # Construct necessary paths and variables
    echo "cellbender count"
    data_dir="${base_data_dir}/${donor}"

    # Find the sample_raw_feature_bc_matrix.h5 file in subdirectories
    input=$(find "$data_dir" -name "raw_feature_bc_matrix.h5" | head -n 1)

    # Check if the file was found
    if [ -z "$input" ]; then
        echo "Error: sample_raw_feature_bc_matrix.h5 not found in any subdirectory of $data_dir"
        exit 1
    fi

    out="${data_dir}/${donor}.cellbender.h5"
    checkpoint="ckpt.tar.gz"

    # Ensure the output directory is writable
    if [ ! -w "$data_dir" ]; then
        echo "Error: Cannot write to directory $data_dir. Ensure the directory is write accessible."
        exit 1
    fi

    # Calculate total droplets included
    if (( expectedCell > 20000 )); then
        total_droplets_included=$(( expectedCell * 13 / 10 ))
    else
        total_droplets_included=20000
    fi

    # Construct cellbender command with the --checkpoint argument
    command=("cellbender" "remove-background"
        "--input" "$input"
        "--output" "$out"
        "--checkpoint" "$checkpoint"
        "--cuda"
        "--expected-cells" "$expectedCell"
        "--total-droplets-included" "$total_droplets_included"
        "--fpr" "0.07"
        "--learning-rate" "1e-5"
    )

    # Print and execute the command
    echo "Running command for donor $donor:"
    echo "${command[@]}"
    "${command[@]}"

    # Check if the checkpoint file was created
    if [ ! -f "$checkpoint" ]; then
        echo "Error: Checkpoint file $checkpoint was not created."
        exit 1
    fi
elif [ "$pipeline_used" == "multi" ]; then
        # Construct necessary paths and variables
    echo "cellbender multi"
    data_dir="${base_data_dir}/${donor}/outs/multi/count"

    # Find the sample_raw_feature_bc_matrix.h5 file in subdirectories
    input=$(find "$data_dir" -name "raw_feature_bc_matrix.h5" | head -n 1)

    # Check if the file was found
    if [ -z "$input" ]; then
        echo "Error: sample_raw_feature_bc_matrix.h5 not found in any subdirectory of $data_dir"
        exit 1
    fi

    out="${data_dir}/${donor}.cellbender.h5"
    checkpoint="ckpt.tar.gz"

    # Ensure the output directory is writable
    if [ ! -w "$data_dir" ]; then
        echo "Error: Cannot write to directory $data_dir. Ensure the directory is write accessible."
        exit 1
    fi

    # Calculate total droplets included
    if (( expectedCell > 20000 )); then
        total_droplets_included=$(( expectedCell * 13 / 10 ))
    else
        total_droplets_included=20000
    fi

    # Construct cellbender command with the --checkpoint argument
    command=("cellbender" "remove-background"
        "--input" "$input"
        "--output" "$out"
        "--checkpoint" "$checkpoint"
        "--cuda"
        "--expected-cells" "$expectedCell"
        "--total-droplets-included" "$total_droplets_included"
        "--fpr" "0.07"
        "--learning-rate" "1e-5"
    )

    # Print and execute the command
    echo "Running command for donor $donor:"
    echo "${command[@]}"
    "${command[@]}"

    # Check if the checkpoint file was created
    if [ ! -f "$checkpoint" ]; then
        echo "Error: Checkpoint file $checkpoint was not created."
        exit 1
    fi
else 
    echo "The pipeline isn't valid"
fi