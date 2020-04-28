The directory structure is:

bin/
The makefile places executable here.  There are aosl example model parameter files, model output and bash shell wrapper scripts here

driver/
The driver code and Makefile are here.

namelist/ 
This directory contains example namelist files.  The file namelist.model.example  Has comments for every line of the input file giving a brief description of what each variable is.

sac/
All the SAC model code is here

share/
Contains routines used by code to read input files, do other things.  Still needs some cleaning for unused code...

snow19/
Snow-17 code is here.  Directory is snow19 because it is an updated version of the original Snow-17, but its still called Snow-17.



The model is run from the command line as:

./lump_model.exe namelist.model

The namelist.model can have any name you want, the code takes this as a command line argument and reads that file as the namelist.
The namelist should contain all the information found in the example namelist for the model to run.  It'll run very quickly and generate output at the file specified by model_out in the namelist.

