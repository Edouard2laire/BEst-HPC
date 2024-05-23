#!/bin/bash 

#module load matlab/R2019b

cd ~/Documents/Project/wMEM-fnirs
matlab -nodisplay -nosplash -nodesktop -r "run_MEM('$1','$2')"
