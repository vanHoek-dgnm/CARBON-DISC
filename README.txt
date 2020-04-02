README to run CARBON-DISC 1.0

COMPUTATIONAL REQUIREMENTS:
LINUX OS
Python 3.5 or more
Anaconda 3 environment (https://www.anaconda.com/distribution/)
Installed numpy module (https://numpy.org/)
Installed NETCDF4 module (https://pypi.org/project/netCDF4/)
Installed mocsy module (http://ocmip5.ipsl.jussieu.fr/mocsy/)

LICENSE:
# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

STEPS TO RUN MODEL:
1. Copy the following directories to an empty directory:
	- A_source_code
	- B_model_input

2. Open terminal 

3. Navigate to directory "core" in directory "A_source_code"

4a. Run abiotic run with command:		
	python dgnm_main.py --inifile ../ini/cmd_m_50yrs_abio_def.ini

4b. Run respiration run with command: 	
	python dgnm_main.py --inifile ../ini/cmd_m_50yrs_resp_def.ini

4c. Run biology run with command: 		
	python dgnm_main.py --inifile ../ini/cmd_m_50yrs_bio_def.ini

When the model runs have been succesful, it is essential to convert the produced raw output to NETCDF files with command:

5.  python ../carbon/code/output_conversion.py outputdir* NETCDF
    * fill in the actual output directory that is being set in the inifile

6.  All output (variables, environmental conditions, fluxes) can be accessed in the produced NETCDF files

7.  If wanted, run aggregate_timeseries.py produces csv tables providing basin aggregated values of almost all used and produced model output.*
    The csv tables are found in the ANALYSIS directory within the output directory 
        python aggregate_timeseries.py --inifile ***
	* Running the output_conversion.py script is necessary before running aggreate_timeseries.py
	** Run the script from the directory /A_source_code/carbon/code
	** fill in the inifile that is used for the model run

8. To conduct a shorter modelrun, make a new inifile with a new output directory and change the starttime and endtime (i.e. 1996 and 2000)

