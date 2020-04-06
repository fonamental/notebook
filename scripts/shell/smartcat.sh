#! /bin/bash

# function smartcat(){                          

if [ $# -lt 1 ]
then
    echo "Concatenate a list of files with returns between"
    echo "Usage: smartcat *.fta > output.txt"
   echo "Note, don't use the same extension for output and input"

 else     
    for FN in "$@" 
		do   
#        echo "File $FN" 
		  	tr "\r" "\n" < "$FN" 
#             cat "$FN" 
			 printf "\n"
         done 
fi
