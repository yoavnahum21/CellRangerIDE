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
    elif pipeline_num == 'flex':
        return my_pipe.multi_flex()
    elif pipeline_num == 'velocyto':
        return my_pipe.velocyto()
    elif pipeline_num == 'QC':
        return my_pipe.qcScores()
    else: 
        print("Your choice is not valid!!!\nPlz choose another pip3line :P ")
    return None


my_pipe = pipeline()

which_pipeline(my_pipe.pipeline)
print("Files were generated successfully")
print("Ciao!")


