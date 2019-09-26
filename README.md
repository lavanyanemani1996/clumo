This is the the description of the repository "clumo". It was created as part of a master thesis project at Argelander Institute for Astronomy, University of Bonn under the supervision of Prof. Dr. Thomas Reiprich, Dr. Florian Pacaud & Dr. Aarti Nagarajan.

clumo can be used to create a mock catalog consisting of various galaxy cluster properties like luminosity, temperature, y-compton paramater, fluxes and x-ray cts. For this the user can select cosmology, sample size, scaling relations and instrument's response. 

1. import cluster_run as cr
2. Chose a name for the directory where you would save your results: For example: name = 'run_1' 
2. Use cr.create_dir(name).
3. Change parameters as required in the parameter_files folder (i.e. clumo/run_1/parameter_files) inside your created directory. NOTE: Do not change parameters inside the default parameter files (i.e. clumo/paramete_files)
4. Use cr.clumorun(name, numberofcores) to create the sample. NOTE: Recommended using numberofcores >= 4

