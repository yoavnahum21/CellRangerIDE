import pandas as pd
import os
import anndata as ad
import torch
import subprocess
import logging
import shutil
import glob

# Configure logging
logging.basicConfig(level=logging.INFO)

# Check if CUDA is available and initialize the CUDA device
def check_cuda_availability():
    if torch.cuda.is_available():
        device = torch.device('cuda')
        logging.info("CUDA is available!")
        logging.info(f"CUDA Device: {torch.cuda.get_device_name(torch.cuda.current_device())}")
        logging.info(f"CUDA Version: {torch.version.cuda}")
    else:
        device = torch.device('cpu')
        logging.info("CUDA is not available. Using CPU instead.")
    return device

def run_cellbender_shell(donor, base_data_dir, script_sh, pipeline_used):
    data_dir = os.path.join(base_data_dir,donor)
    if pipeline_used == 'count':
        # Construct file path to metrics summary CSV
        
        file_pattern = os.path.join(data_dir, "metrics_summary.csv")
        
        # Read CSV file and extract expectedCell value
        try:
            df = pd.read_csv(file_pattern)
            print("DataFrame read from CSV:")
            print(df)
            # Extract the 'Metric Value' where 'Metric Name' is 'Cells'
            expectedCell = df.loc[df['Metric Name'] == 'Cells', 'Metric Value'].values[0]
            print(f"Extracted expectedCell value: {expectedCell}")
            # Remove commas from the expectedCell value
            expectedCell = int(expectedCell.replace(',', ''))
            print(f"Processed expectedCell value: {expectedCell}")
        except Exception as e:
            print(f"Error reading expectedCell from CSV: {e}")
            return False
        
        script_path = script_sh
        
        # Command to execute the shell script
        command = [
            "bash",
            script_path,
            donor,
            base_data_dir,
            str(expectedCell),
            pipeline_used  
        ]

        print(f"Running shell script for donor {donor}...")
        try:
            subprocess.run(command, check=True)
            print(f"Successfully ran shell script for donor {donor}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error running shell script for donor {donor}: {e}")
            return False
        
    if pipeline_used == 'multi':
        expectedCell = 0
        file_pattern = os.path.join(data_dir, "outs/per_sample_outs")
        for sample_dir in os.listdir(file_pattern):
            if sample_dir.startswith('.'):
                continue
            sample_pattern = os.path.join(file_pattern, sample_dir)
            metrics_path = os.path.join(sample_pattern, "metrics_summary.csv")
            
            try:        
                df = pd.read_csv(metrics_path)
                print("DataFrame read from CSV:")
                print(df)
                # Extract the 'Metric Value' where 'Metric Name' is 'Cells'
                countcells = df.loc[df['Metric Name'] == 'Cells', 'Metric Value'].values[0]
                print(f"Extracted expectedCell value: {expectedCell}")
                # Remove commas from the expectedCell value
                expectedCell += int(countcells.replace(',', ''))
                print(f"Processed expectedCell value: {expectedCell}")
            
            except Exception as e:
                print(f"Error reading expectedCell from CSV: {e}")
                return False
            
        script_path = script_sh
        
        # Command to execute the shell script
        command = [
            "bash",
            script_path,
            donor,
            base_data_dir,
            str(expectedCell),
            pipeline_used  
        ]

        print(f"Running shell script for donor {donor}...")
        try:
            subprocess.run(command, check=True)
            print(f"Successfully ran shell script for donor {donor}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error running shell script for donor {donor}: {e}")
            return False

            


def save_files_locally(donor, base_data_dir, save_dir):
    data_dir = os.path.join(base_data_dir, donor)

    # Use glob to find the raw_feature_bc_matrix.h5 file in all subdirectories
    pattern = os.path.join(data_dir, "**", "raw_feature_bc_matrix.h5")
    file_list = glob.glob(pattern, recursive=True)

    # Check if the file was found
    if file_list:
        input = file_list[0]  # Take the first match if there are multiple
        print(f"Found file: {input}")
    else:
        print("raw_feature_bc_matrix.h5 not found in any subdirectory")
        return False
    
    # Ensure the save directory exists and is writable
    if not os.path.exists(save_dir):
        try:
            os.makedirs(save_dir)
            print(f"Created save directory: {save_dir}")
        except Exception as e:
            print(f"Error creating save directory {save_dir}: {e}")
            return False
    elif not os.access(save_dir, os.W_OK):
        print(f"Cannot write to save directory {save_dir}. Ensure the directory is write accessible.")
        return False
    
    # Define the paths for the output files
    out = os.path.join(data_dir, f"{donor}.cellbender.h5")
    out5H = os.path.join(data_dir, f"{donor}.cellbender.h5ad")
    # ckpt_file = os.path.join(data_dir, f"{donor}.ckpt.tar.gz")

    # Debugging: Print the paths
    print(f"Output file: {out}")
    # print(f"Checkpoint file: {ckpt_file}")
    
    try:
        # Load the input H5 file directly into an AnnData object
        adata = ad.read_h5ad(input)
        adata.uns.clear()

        # Save the processed AnnData object back to a new H5AD file
        adata.write_h5ad(out5H)
        print(f"Processed donor {donor}. Output saved to {out5H}")
    except Exception as e:
        print(f"Failed to process output for donor {donor}: {str(e)}")
        return False

    files_to_save = [
        f"{donor}.ckpt.tar.gz",
        f"{donor}.cellbender.h5ad",
        f"{donor}.cellbender_filtered.h5",
        f"{donor}.cellbender_report.html",
        f"{donor}.cellbender.log",
        f"{donor}.cellbender.pdf"
    ]

    for file_name in files_to_save:
        file_path = os.path.join(data_dir, file_name)

        # Debugging: Print the file path being checked
        print(f"Checking file: {file_path}")
        if os.path.exists(file_path):
            try:
                dest_path = os.path.join(save_dir, file_name)
                shutil.copy(file_path, dest_path)
                print(f"Saved {file_name} to {dest_path}")
            except Exception as e:
                print(f"Failed to save {file_name}. Error: {e}")
        else:
            print(f"File {file_name} does not exist. Skipping save.")

    return True

