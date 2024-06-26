#!/bin/bash

Help()
{
   # Display Help
   echo "Invalid Usage"
   echo "Syntax: launch_MEM.sh token OPTIONS.json"
   echo "options:"
   echo "token            Token of the data you want to process."
   echo "OPTIONS.json     Options used for MEM."
   echo
}

if (( "$#" != 2 )) 
then
    Help
    exit;
fi


token=${1}
option=${2}

mkdir data/${token}/logs 

iFile=1
for FILE in data/${token}/in/*
do
	echo $FILE;
	qsub -j y -o data/${token}/logs/${iFile}_${token}.txt \ 
		 -pe smp 15 -S /bin/bash  -cwd -q  all.q -N MEM_${iFile}_${token} \
		./run_MEM_HPC.sh $FILE ${option}
    iFile=${iFile}+1
done
