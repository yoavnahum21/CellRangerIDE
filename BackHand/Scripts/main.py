from pipeline import pipeline
from _utils import PROGRAM_PATH
import pandas as pd
import os 


def which_pipeline(pipeline_num) -> None:
    
    if pipeline_num == '1':
        my_pipe.mkfastq()    
    elif pipeline_num == '2':
        my_pipe.multiplex()
    elif pipeline_num == '3':
        my_pipe.count()
    elif pipeline_num == '4':
        my_pipe.demultiplex()
    elif pipeline_num == '5':
        my_pipe.run_basic_pipeline()
    elif pipeline_num == '6':
        return None
    else: 
        print("Your choise is not valid!!!\nPlz choose another one ")
    return None


# csv_address = ("/home/labs/nyosef/yoavnah/CellRangerIDE/Projects/yoav/File_Path/CSV/program.csv")
# csv_file = pd.read_csv(csv_address)

my_pipe = pipeline()
print("Which pipeline would you like to use??\n 1) Mkfastq\n 2) Multiplexing\n 3) Count\n 4) Demultiplexing\n 5) Run the wholllllleeee pipeline with a basic example\n 6) End task")
pipe = input()
which_pipeline(pipe)
while pipe != '6':
    print("Which pipeline would you like to use now??\n 1) Mkfastq\n 2) Multiplexing\n 3) Count\n 4) Demultiplexing\n 5) Run the wholllllleeee pipeline with a basic example\n 6) End task")
    pipe = input()
    which_pipeline(pipe)
    
print("Ciao")
