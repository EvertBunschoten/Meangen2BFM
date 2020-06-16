
import numpy as np


class StagenReader:
    def __init__(self):
        row = 1
        sec = 1



        section = False
        with open("stagen.dat", "r") as stagen:

            for line in stagen:
                s = " ".join(line.split())
                z = s.split(" ")
                if len(z) == 5 and z[-2] + z[-1] == 'NSECTIONS':
                    nRows = int(z[0])
                    nSec = int(z[1])

                    theta_in = np.zeros([nSec, nRows])
                    theta_out = np.zeros([nSec, nRows])
                    t_le = np.zeros([nSec, nRows])
                    t_te = np.zeros([nSec, nRows])
                    t_max = np.zeros([nSec, nRows])
                    x_tmax = np.zeros([nSec, nRows])
                    X_LE = np.zeros([nSec, nRows])
                    Z_LE = np.zeros([nSec, nRows])
                    X_TE = np.zeros([nSec, nRows])
                    Z_TE = np.zeros([nSec, nRows])

                    section = True
                if section:
                    if sec <= nSec:
                        if len(z) == 6:
                            if z[0] == '0.0000' and z[-4] + z[-3] + z[-2] + z[-1] == 'BLADECENTRELINEANGLES':
                                theta_in[sec-1, row-1] += float(z[1])
                            elif z[0] == '1.0000' and z[-4] + z[-3] + z[-2] + z[-1] == 'BLADECENTRELINEANGLES':
                                theta_out[sec-1, row-1] += float(z[1])
                               # sec += 1
                        if len(z) == 10 and z[-3] + z[-2] + z[-1] == 'BLADEPROFILESPECIFICATION':
                            t_max[sec-1, row-1] += float(z[2])
                            x_tmax[sec-1, row-1] += float(z[3])
                            t_le[sec-1, row-1] += float(z[4])
                            t_te[sec-1, row-1] += float(z[5])
                        if len(z) == 9 and z[-2] + z[-1] == 'EDGECOORDINATES':
                            X_LE[sec-1, row-1] += float(z[0])
                            Z_LE[sec-1, row-1] += float(z[2])
                            X_TE[sec - 1, row - 1] += float(z[1])
                            Z_TE[sec - 1, row - 1] += float(z[3])
                            sec += 1

                    if sec == nSec + 1:
                        sec = 1
                        row += 1
        self.theta_in = theta_in
        self.theta_out = theta_out
        self.t_max = t_max
        self.x_tmax = x_tmax
        self.t_le = t_le
        self.t_te = t_te
        self.X_LE = X_LE
        self.X_TE = X_TE
        self.Z_LE = Z_LE
        self.Z_TE = Z_TE