#!/home/evert/anaconda3/bin/python3

# --------------------------------------------------------------------------------------------------------------- #
# This code serves the purpose of using general engineering inputs for axial turbomachinery design and translate them
# to a suitable 2 or 3-dimensional mesh and body-force method input for SU2. This program works in conjunction with
# Meangen, Parablade(the body-force branch), which generates the detailed blade shape from Meangen output.
# For 3D meshing, Gmesh has to be installed and added to the PYTHONPATH. For 2D meshing, UMG2 has to be installed and
# added to the PATH as well. For body-force analyses in SU2, the feature_bodyforce_turbo has to be cloned and installed
# on the user's machine. Template files for all relevant configuration files can be found in the 'templates' folder and
# all scripts are located in the 'executables' folder.
#
# Author: E.C.Bunschoten
# Institution: TU Delft
# Date: 15/6/2020
# --------------------------------------------------------------------------------------------------------------- #
import sys
import os
import time

# Getting the executables directory
HOME = os.environ["M2BFM"]
sys.path.append(HOME + "executables/")
# Importing all relevant executables from the installation directory
from Meangen2Parablade import Meangen2Parablade
from Parablade2UMG2 import WriteUMG, writeStageMesh_BFM, writeStageMesh_Blade
from SU2Writer import writeBFMinput, ReadUserInput, writeSU2input
from Mesh3D import Gmesh3D, Gmesh2D, FullAnnulus
from dataPlotter import axial_data_plotter
# from ParaviewPost import AxialMachine

# Reading input file
DIR = os.getcwd() + '/'
try:
    INFile = DIR + sys.argv[-1]
except:
    INFile = DIR + 'M2P.cfg'      # Default File name
try:
    IN = ReadUserInput(INFile)
except:
    raise Exception('\n\n\n''Something went wrong when reading the configuration file,exiting the program...'
                    '\n\nTo call MakeBlade.py from terminal type:'
                    '\n\tMakeBlade.py <configuration file name>')
t_start = time.time()


# Executing Meangen and writing Parablade input files.
M = Meangen2Parablade(IN)

# Extracting stage count and calculating row count.
n_stage = int(IN["N_stage"][0])
n_rows = 2*n_stage

# Looping over the number of stages to create folders for each stage and respective bladerow.
for i in range(n_stage):
    if os.path.isdir("Stage_"+str(i+1)):
        if os.path.isdir("Stage_"+str(i+1) + "/Bladerow_1"):
            os.system("mv Bladerow_"+str(2*i + 1) + ".cfg" + " Stage_"+str(i+1)+"/Bladerow_1/Bladerow.cfg")
        else:
            os.system("mkdir "+"Stage_"+str(i+1)+"/Bladerow_1")
            os.system("mv Bladerow_" + str(2 * i + 1) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_1/Bladerow.cfg")
        if os.path.isdir("Stage_"+str(i+1) + "/Bladerow_2"):
            os.system("mv Bladerow_" + str(2 * i + 2) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_2/Bladerow.cfg")
        else:
            os.system("mkdir " + "Stage_" + str(i + 1) + "/Bladerow_2")
            os.system("mv Bladerow_" + str(2 * i + 2) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_2/Bladerow.cfg")
    else:
        os.system("mkdir Stage_"+str(i+1))
        os.system("mkdir " + "Stage_" + str(i + 1) + "/Bladerow_1")
        os.system("mv Bladerow_" + str(2 * i + 1) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_1/Bladerow.cfg")
        os.system("mkdir " + "Stage_" + str(i + 1) + "/Bladerow_2")
        os.system("mv Bladerow_" + str(2 * i + 2) + ".cfg" + " Stage_" + str(i + 1) + "/Bladerow_2/Bladerow.cfg")

# Checking for body-force and/or blade mesh option.
BFM = False
Blade = False
if IN["MESH_BFM"] == 'YES':
    BFM = True
if IN["MESH_BLADE"] == 'YES':
    Blade = True

# Looping over the blade rows to execute Parablade for each blade row and writing UMG2 input files for each blade row
# mesh. These individual meshes will be combined later into a complete machine mesh.
row = 1
for i in range(n_stage):
    for j in [1, 2]:
        # Moving to current blade row directory.
        os.chdir(os.getcwd()+"/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/")

        # Checking for blade plot option.
        if IN['PLOT_BLADE'] == 'YES':
            print("plotting blade")
            os.system("PlotBlade.py Bladerow.cfg > Parablade.out")

        # Executing Parablade.
        print("Running Parablade...")
        os.system("MakeBlade.py Bladerow.cfg > Parablade.out")
        print("Done!")

        if IN["ADJOINT"] == 'YES':
            os.system("sed -i 's/GEOMETRY/SENSITIVITY/g' Bladerow.cfg")
            print("Getting sensitivities for blade row "+str(2*i + j))
            os.system("MakeBlade.py Bladerow.cfg")
            print("Done!")
        # In case the dimension number is 2, UMG2 input files will be written, depending on the mesh case option.
        #if IN['N_dim'][0] == 2:
            #WriteUMG(j, i+1, M, IN, bodyForce=BFM, blade=Blade)

        # Updating row count.
        row += 1
        os.chdir(DIR)

# The individual 2D meshes are combined into a full machine mesh. In case the dimension number is 3 and a BFM mesh is
# desired, a suitable 3D mesh will be written. This option is currently not yet available for physical blades.

if BFM:
    # Writing BFM input file suitable for SU2 BFM analysis.
    print("Writing Body-force SU2 input file...", end='     ')
    writeBFMinput(M)
    print("Done!")

    # Writing 3D BFM mesh or combining individual 2D blade row meshes depending on case dimension.
    if IN['N_dim'][0] == 3:
        print("Writing 3D BFM mesh:...")
        Gmesh3D(M, IN)
        print("Done!")
    else:
        print("Writing 2D BFM mesh...", end='     ')
        Gmesh2D(M, IN)
        print("Done!")

#
if Blade:
    if IN['N_dim'][0] == 3:
        print("3D physical blade meshing is not yet implemented!")
    else:
        os.chdir(DIR)
        print("Writing 2D Blade analysis SU2 machine mesh file...", end='     ')
        writeStageMesh_Blade(M)
        print("Done!")
print("Total geometry and mesh generation took "+str(format(time.time() - t_start, ".2f")) + " seconds")
writeSU2input(IN)

if IN["SOLVE"] == 'AUTOMATIC':
    os.system("SU2_CFD BFM_comp.cfg")

if IN["POSTPROCESS"] == 'YES':
    print("Postprocessing simulation data....")
    os.system("pvpython "+HOME+"executables/ParaviewPost.py "+INFile)
    print("done!")
    print("Creating data plots...")
    axial_data_plotter()
    print("Done!")