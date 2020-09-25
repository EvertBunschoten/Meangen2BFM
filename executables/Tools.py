import numpy as np
from scipy.interpolate import griddata

class RaycastInterpolation:
    def __init__(self, inFile, x, r):
        self.f = open(inFile, "r")
        self.X_field = x
        self.R_field = r
        infoLine = self.f.readline().strip().split('\t')
        self.n_blades = int(infoLine[0])
        self.n_points = int(infoLine[1])
        self.n_sec = int(infoLine[2])
        self.inFile = inFile
        self.GetPrimaryMatrix()
        self.GetGriddata()
        self.f.close()

    def GetGriddata(self):
        self.Interpvalues = np.zeros((len(self.X_field), 3))
        for i in range(self.n_blades):
            X = self.x_primary[i].flatten()
            R = self.r_primary[i].flatten()
            points = np.array([X, R]).transpose()
            values_dNx_dphi = self.Nx_primary[i].flatten()
            values_dNt_dphi = self.Nt_primary[i].flatten()
            values_dNr_dphi = self.Nr_primary[i].flatten()
            dNx_interp = griddata(points=points, values=values_dNx_dphi, xi=(self.X_field, self.R_field), method='linear', fill_value=0.0)
            dNt_interp = griddata(points=points, values=values_dNt_dphi, xi=(self.X_field, self.R_field), method='linear', fill_value=0.0)
            dNr_interp = griddata(points=points, values=values_dNr_dphi, xi=(self.X_field, self.R_field), method='linear', fill_value=0.0)

            self.Interpvalues[:, 0] += np.transpose(dNx_interp)
            self.Interpvalues[:, 1] += np.transpose(dNt_interp)
            self.Interpvalues[:, 2] += np.transpose(dNr_interp)

    def GetPrimaryMatrix(self):

        self.x_primary = np.zeros((self.n_blades, self.n_sec, self.n_points))
        self.r_primary = np.zeros((self.n_blades, self.n_sec, self.n_points))
        self.Nx_primary = np.zeros((self.n_blades, self.n_sec, self.n_points))
        self.Nt_primary = np.zeros((self.n_blades, self.n_sec, self.n_points))
        self.Nr_primary = np.zeros((self.n_blades, self.n_sec, self.n_points))

        f = open(self.inFile, "r")
        lines = f.readlines()
        s = 1
        self.x_min = 0.0
        for i in range(self.n_blades):
            for j in range(self.n_sec):
                s += 1
                for k in range(self.n_points):
                    line = lines[s].strip().split('\t')
                    self.x_primary[i, j, k] = float(line[0])
                    self.r_primary[i, j, k] = float(line[1])
                    self.Nx_primary[i, j, k] = float(line[2])
                    self.Nt_primary[i, j, k] = float(line[3])
                    self.Nr_primary[i, j, k] = float(line[4])
                    if float(line[0]) < self.x_min:
                        self.x_min = float(line[0])
                    s += 1

        f.close()