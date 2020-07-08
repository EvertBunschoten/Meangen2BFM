import numpy as np
import os
import sys
# Getting the executables directory
HOME = os.environ["M2BFM"]
sys.path.append(HOME + "executables/")

# Importing all relevant executables from the installation directory
from Meangen2Parablade import Meangen2Parablade
class dN_dalpha:
    def __init__(self, IN):
        self.Variables = IN["DESIGN_VARIABLES"]

        self.step = 1e-3
        self.IN_new = IN
        self.n_rows = int(IN["N_stage"][0]) * 2
        self.Mesh_points = int(IN["AXIAL_POINTS"][0])
        os.system("mkdir ./PartialGradients")
        self.CentralDifferences()
    def CentralDifferences(self):
        if isinstance(self.Variables, str):
            L = 1
        else:
            L = len(self.Variables)

        for i in range(L):
            if isinstance(self.Variables, str):
                var = self.Variables
            else:
                var = self.Variables[i]

            print("Computing camber normal derivative of design variable: " + var)
            for j in range(len(self.IN_new[var])):

                for k in range(self.n_rows):
                    self.IN_new[var][j] += self.step
                    M = Meangen2Parablade(self.IN_new)

                    Normals = self.GetCamberNormals(k+1)
                    N_plus = Normals[:, 2:]
                    XR_plus = Normals[:, :2]

                    self.IN_new[var][j] -= 2 * self.step
                    M = Meangen2Parablade(self.IN_new)

                    Normals = self.GetCamberNormals(k+1)
                    N_minus = Normals[:, 2:]
                    XR_minus = Normals[:, :2]

                    if k == 0:
                        self.dN_dx = np.concatenate((0.5 * (XR_plus + XR_minus), (N_plus - N_minus)/(2 * self.step)), axis=1)
                    else:
                        self.dN_dx = np.append(self.dN_dx, np.concatenate((0.5 * (XR_plus + XR_minus), (N_plus - N_minus)/(2 * self.step)), axis=1), axis=0)

                varName = var

                self.WriteDerivatives(dN_dx=self.dN_dx, varName=varName)
                self.IN_new[var][j] += self.step

    def GetCamberNormals(self, rowNumber):
        os.system("MakeBlade.py Bladerow_"+str(rowNumber)+".cfg > thingy")
        BFMfile = open("./output/mesh_files/BFM_input", "r")

        N = np.zeros((self.Mesh_points**2, 8))
        j = 1
        s = 0
        Lines = BFMfile.readlines()
        for i in range(self.Mesh_points):
            j += 1
            #print(j)
            lines = Lines[j:j + self.Mesh_points]
            #print(lines)
            for k in range(len(lines)):
                line = lines[k].strip().split('\t')
                #print(line)

                N[s, :] += np.array([float(l) for l in line])
                #print(N[s, -1])
                s += 1

            j += self.Mesh_points
        BFMfile.close()
        #os.system("rm -rf ./output")
        return N

    def WriteDerivatives(self, dN_dx, varName):
        file = open("./PartialGradients/"+varName, "w+")
        file.write(str(self.n_rows)+"\t"+str(self.Mesh_points)+"\t"+str(self.Mesh_points)+"\n")
        #file.write("# x, r, dNx_dvar, dNt_dvar, dNr_dvar\n")

        sec_count = 0
        for i in range(self.n_rows*self.Mesh_points**2):
            if i % self.Mesh_points == 0:
                file.write("Section " + str(sec_count) + "\n")
                file.write("\t".join([str(s) for s in dN_dx[i, :]])+"\n")
                sec_count += 1
            else:
                file.write("\t".join([str(s) for s in dN_dx[i, :]]) + "\n")
        file.close()




