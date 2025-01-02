# 
from pipeline import pipeline

def which_pipeline(pipeline_num) -> None:
    
    if pipeline_num == 'mkfastq':
        my_pipe.mkfastq()    
    elif pipeline_num == 'multi':
        my_pipe.multiplex()
    elif pipeline_num == 'count':
        my_pipe.count()
    elif pipeline_num == 'demulti':
        my_pipe.demultiplex()
    elif pipeline_num == 'mkref':
        my_pipe.make_custom_reference()
    elif pipeline_num == 'cellbender':
        my_pipe.cellbender()
    elif pipeline_num == '7':
        return None
    else: 
        print("Your choise is not valid!!!\nPlz choose another pip3line :P ")
    return None


my_pipe = pipeline()

which_pipeline(my_pipe.pipeline)
print("Ciao!")
