#!/home/evert/anaconda3/bin/python3
import sys
import os


import time

HOME = os.environ["M2BFM"]
sys.path.append(HOME + "executables/")

from Meangen2Parablade import Meangen2Parablade
from Parablade2UMG2 import WriteUMG, writeStageMesh, writeStageMesh_Blade
from SU2Writer import writeBFMinput, ReadUserInput, writeSU2input
from Mesh3D import Gmesh3D

DIR = os.getcwd() + '/'


try:
    IN = ReadUserInput(DIR + 'M2P.cfg')
except:
    raise Exception('\n\n\n''Something went wrong when reading the configuration file,exiting the program...'
                    '\n\nTo call MakeBlade.py from terminal type:'
                    '\n\tMakeBlade.py <configuration file name>')
t_start = time.time()

M = Meangen2Parablade(IN)

n_stage = int(IN["N_stage"][0])
n_rows = 2*n_stage


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

row = 1
BFM = True
Blade = False


for i in range(n_stage):
    for j in [1, 2]:
        os.chdir(os.getcwd()+"/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/")
        if IN['PLOT_BLADE'] == 'YES':
            print("plotting blade")
            os.system("PlotBlade.py Bladerow.cfg")
        os.system("MakeBlade.py Bladerow.cfg")

        if IN['N_dim'][0] == 2:
            WriteUMG(j, i+1, M, bodyForce=BFM, blade=Blade)
            # print("Starting UMG2 meshing...", end='                 ')
            # os.system("HYMESH.sh > mesher.log")
            # print("Done!")

        row += 1
        os.chdir(DIR)
    # os.chdir(DIR+"Stage_"+str(i))

if BFM:
    print("Writing Body-force SU2 input file...", end='     ')
    writeBFMinput(M)
    print("Done!")

    if IN['N_dim'][0] == 3:
        print("Writing 3D BFM mesh:...")
        Gmesh3D(M)
        print("Done!")
    else:
        print("Writing Body-force SU2 machine mesh file...", end='     ')
        writeStageMesh(M)
        print("Done!")
# if Blade:
#     os.chdir(DIR)
#     print("Writing Blade analysis SU2 machine mesh file...", end='     ')
#     writeStageMesh_Blade(M)
#     print("Done!")
# print("Total geometry and mesh generation took "+str(format(time.time() - t_start, ".2f")) + " seconds")
# writeSU2input()
