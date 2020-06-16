Meangen2BFM is a program which allows for the creation of 2D and 3D body-force and (future) blade computation meshes for axial turbomachinery analysis from Meangen input. The user specifies the machine type, number of stages, duty coefficients, performance parameters and blade thickness parameters and the program creates a 2D or 3D mesh, alongside the input file required for body-force analysis in SU2. 

For the program to work, installation of Parablade is required(the 'body-force' branch), GMesh and UMG2 and their installation location has to have been added to PATH and PYTHONPATH. There is a copy of Meangen already present in the 'executables' folder. The following lines should be added to the .bashrc file:

export M2BFM="<path to Meangen2BFM folder>"
export PATH=$PATH:"<path to Meangen2BFM folder>/executables/"
export PYTHONPATH=$PYTHONPATH:"<path to Meangen2BFM folder>/executables/"

A template of the inputfile to be used can be found in the 'templates' folder under the name 'M2P.cfg'. More detailed explanations on all the input options can be found in there. There is also a separate template file containing the inputs for the meshing part of the program. This file can also be found in the 'templates' folder under the name 'Mesh_params.cfg'. 

To execute the program, the user has to have the  file and the 'Mesh_params.cfg' file in the target directory. By then executing the following command in the terminal,

Meangen2BFM.py inputfile.cfg

the program will initiate and start generating output folders for Meangen and Parablade, alongside a mesh suitable for body-force analysis or blade computation. Currently, full integration with SU2 is not yet functional, so SU2 has to be initiated manually. A template of the SU2 configuration file required for the body-force analysis can be found in the 'templates' file under the name 'BFM_comp.cfg'. 


