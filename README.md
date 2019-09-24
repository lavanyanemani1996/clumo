This is the the description of the repository "clumo". It was created as part of a master thesis project at Argelander Institute for Astronomy, University of Bonn under the supervision of Prof. Dr. Thomas Reiprich in 2019. 

clumo can be used to create a mock catalog consisting of various galaxy cluster properties like luminosity, temperature, y-compton paramater, fluxes and x-ray cts. For this the user can select cosmology, sample size, scaling relations and instrument's response. 

1. import cluster_run as cr
2. Use cr.create_dir() to create a directory where you would save your results.
3. Change parameters as required in the parameter_files folder inside your created directory
4. Use cr.clumorun(directoryname, numberofcores) to create the sample
