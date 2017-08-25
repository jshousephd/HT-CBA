#############################################################################################################
##Author :          John House - BRC;NCSU
##Last Modified :   2/28/17 - John House
##Purpose:          Create Directory Structures For New Project 
##
##                  
#############################################################################################################
### Change following path to the base directory of your TempOseq experiments
mainDir = "C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline"
### change following to new project name and execute rest of script ###
subDir = "ht_cba_manuscript"
setwd(mainDir); print(getwd())

dir.create(file.path(mainDir,subDir))
dir.create(file.path(mainDir,subDir,"Input Files"))    ## Sequencing Files and Hash File for Experiment
dir.create(file.path(mainDir,subDir,"R Scripts"))      ## All other production Scripts
dir.create(file.path(mainDir,subDir,"Output Files"))   ## 
dir.create(file.path(mainDir,subDir,"Output Files/MaxDoseContrasts"))
dir.create(file.path(mainDir,subDir,"Output Files/DRM"))
dir.create(file.path(mainDir,subDir,"Output Files/Cleaned_Data"))
dir.create(file.path(mainDir,subDir,"Output Files/RawCountResults"))
dir.create(file.path(mainDir,subDir,"Output Files/DESeq2 dds Files"))


