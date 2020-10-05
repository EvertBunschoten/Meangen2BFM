# --------------------------------------------------------------------------------------------------------------- #
# This program writes the input file for the BFM solver in SU2. It does so by combining the BFM input files
# written for each blade row by Parablade. In the future, writeSU2input class will be able to write the input file for
# SU2.
# --------------------------------------------------------------------------------------------------------------- #

import os
import re
import numpy as np

# Class used for combining all individual blade row BFM input files into a single file for the entire machine.
class writeBFMinput:
    n_sec =     None  # Number of sections in spanwise direction.
    n_points =  None  # Number of points in axial direction.
    n_rows =    0     # Number of blade rows.

    def __init__(self, Meangen):
        # Storing Meangen class input.
        self.M = Meangen

        # Creating machine input file for BFM interpolation.
        self.dir = os.getcwd()
        self.BFMfile = open(self.dir + "/BFM_stage_input", "w+")
        # Getting section and axial point information from Parablade output.
        self.getBFMinputs()
        # Writing number of rows, sections and axial points into the BFM interpolation file.
        self.BFMfile.write("%i\t%i\t%i\n" % (self.n_rows, self.n_sec, self.n_points))
        # Appending the individual blade row BFM files into the main BFM interpolation file.
        self.writeInputFile()
        self.BFMfile.close()

    def getBFMinputs(self):
        # Looping over the number of stages.
        for i in range(self.M.n_stage):
            # Looping over the two blade rows in each stage.
            for j in [1, 2]:
                # Opening the BFM input file of the current blade row.
                with open(self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/output/mesh_files/BFM_input", "r") as file:
                    # Reading axial point and section count information from the first line.
                    first_line = file.readline().split('\t', 2)
                    first_line[-1] = first_line[-1].strip()
                    # Updating axial point count, section count and row count.
                    self.n_sec = int(first_line[0])
                    self.n_points = int(first_line[1])
                    self.n_rows += 1
                file.close()

    def writeInputFile(self):
        # Looping over the number of stages.
        for i in range(self.M.n_stage):
            # Looping over the two blade rows in each stage.
            for j in [1, 2]:
                with open(self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/output/mesh_files/BFM_input", "r") as file:
                    # Skipping the first line in the blade row file.
                    lines = file.readlines()[1:]
                    # Distinguishing between rotor and stator row, depending on machine type.
                    if self.M.machineType == 'C':
                        # For a compressor, the rotor row is followed by a stator row.
                        rotFac = [1, 0]
                        blades = [self.M.N_b_R[i], self.M.N_b_S[i]]
                    else:
                        # For a turbine, the stator row is followed by a rotor row.
                        rotFac = [0, 1]
                        blades = [self.M.N_b_S[i], self.M.N_b_R[i]]
                    for line in lines:
                        # Writing all lines of the blade row BFM file to the machine BFM file, as well as the rotation
                        # factor and the number of blades.
                        self.BFMfile.write(line.strip()+"\t"+str(rotFac[j-1])+"\t"+str(int(blades[j-1]))+"\n")
                file.close()

# Function used to read values from the user input file. This function was copied from Parablade.
def ReadUserInput(name):
    IN = {}
    infile = open(name, 'r')
    for line in infile:
      words = re.split('=| |\n|,|[|]', line)
      if not any(words[0] in s for s in ['\n', '%', ' ', '#']):
        words = list(filter(None, words))
        for i in range(0, len(words)):
            try:
                words[i] = float(words[i])
            except:
                words[i] = words[i]
        if len(words[1::1]) == 1 and isinstance(words[1], float):
            IN[words[0]] = [words[1]]
        elif len(words[1::1]) == 1 and isinstance(words[1], str):
            IN[words[0]] = words[1]
        else:
            IN[words[0]] = words[1::1]
    IN['Config_Path'] = name
    return IN

# This class is still under development. In the future, this class will be able to read the outlet pressure from
# Meangen output and write the SU2 configuration file for BFM or blade analysis.
class writeSU2input:
    fileName = "meandesign.out"
    def __init__(self, IN):
        HOME = os.environ["M2BFM"]
        self.IN = IN
        template_dir = HOME + "templates/"
        os.system("cp "+template_dir+"BFM_comp.template ./BFM_comp.cfg")
        self.ReplaceTerms()

    def ReplaceTerms(self):
        gamma = self.IN["gamma"][0]
        R = self.IN["R_gas"][0]
        rot_axis = self.IN["Rotation_axis"]

        os.system("sed -i 's/GAMMA_FLUID/" + str(gamma) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/R_FLUID/" + str(R) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/CP_FLUID/" + str(gamma * R / (gamma - 1)) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/OMEGA/" + str(self.IN["Omega"][0]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/ROT_X/" + str(rot_axis[0]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/ROT_Y/" + str(rot_axis[1]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/ROT_Z/" + str(rot_axis[2]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/UIN_X/" + str(rot_axis[0]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/UIN_Y/" + str(rot_axis[1]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/UIN_Z/" + str(rot_axis[2]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/MEAN_RADIUS/" + str(self.IN["r_m"][0]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/P_TOT_IN/" + str(self.IN["P_t_in"][0] * 1e5) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/T_TOT_IN/" + str(self.IN["T_t_in"][0]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/P_STAT_OUT/" + str(self.IN["P_s_out"][0] * 1e5) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/WEDGE_X/" + str(rot_axis[0] * self.IN["WEDGE"][0]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/WEDGE_Y/" + str(rot_axis[1] * self.IN["WEDGE"][0]) + "/g' BFM_comp.cfg")
        os.system("sed -i 's/WEDGE_Z/" + str(rot_axis[2] * self.IN["WEDGE"][0]) + "/g' BFM_comp.cfg")
        if self.IN["N_dim"][0] == 3.0:
            os.system("sed -i 's/meshfilename/3DBFM_perio.su2/g' BFM_comp.cfg")
        else:
            os.system("sed -i 's/meshfilename/2DBFM_perio.su2/g' BFM_comp.cfg")
