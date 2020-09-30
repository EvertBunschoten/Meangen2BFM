from paraview.simple import *
import numpy as np

T_in = 281.1614467
P_in = 101478.5239

BFM = LegacyVTKReader(FileNames=["BFMOUTPUT_3.vtk"])

TdS = Calculator(Input=BFM)
TdS.ResultArrayName='dS'
TdS.Function = '(1004.703*ln(Temperature/'+str(T_in)+') - 287.15*ln(Pressure/'+str(P_in)+'))'
RenameSource("dS", TdS)


slice1 = Slice(Input = TdS)
slice1.SliceType.Normal = [0, 0, 1]
slice1.SliceType.Origin = [0, 0, -0.05]

IV = IntegrateVariables(Input = slice1)

zrange = np.linspace(-0.09, 0.15, 1000)

outputFile = open("outputData.txt", "w+")
outputFile.write("ax\tf_p\tTds_dax\n")
ds = []
P = []
BF_z = []
T = []
rho = []
u_x = []
for z in zrange:
    slice1.SliceType.Origin = [0, 0, z]
    iv_data = paraview.servermanager.Fetch(IV)
    area = iv_data.GetCellData().GetArray("Area").GetValue(0)
    BF_z.append(iv_data.GetPointData().GetArray("Body-Force").GetValue(2)/area)
    rho.append(iv_data.GetPointData().GetArray("Density").GetValue(0)/area)
    u_x.append(iv_data.GetPointData().GetArray("Momentum").GetValue(2)/iv_data.GetPointData().GetArray("Density").GetValue(0))
    ds.append(iv_data.GetPointData().GetArray("dS").GetValue(0)/area)
    T.append(iv_data.GetPointData().GetArray("Temperature").GetValue(0)/area)
    P.append(iv_data.GetPointData().GetArray("Pressure").GetValue(0)/area)
    #outputFile.write(str(z)+"\t"+str(BF_z)+"\t"+str(tds)+"\n")
  

Tdsdz = []
dpdz = []
for i in range(len(zrange)):
    if i ==0:
        deltaZ = zrange[i+1] - zrange[i]
        deltaS = ds[i+1] - ds[i]
        deltaP = P[i+1] - P[i]
        deltaMom = rho[i+1]*u_x[i+1]**2 -  rho[i]*u_x[i]**2
    elif i == len(zrange)-1:
        deltaZ = zrange[i] - zrange[i-1]
        deltaS = ds[i] - ds[i-1]
        deltaP = P[i] - P[i-1]
        deltaMom = rho[i]*u_x[i]**2 -  rho[i-1]*u_x[i-1]**2
    else:
        deltaZ = zrange[i+1] - zrange[i-1]
        deltaS = ds[i+1] - ds[i-1]
        deltaP = P[i+1] - P[i-1]
        deltaMom = rho[i+1]*u_x[i+1]**2 -  rho[i-1]*u_x[i-1]**2
     
    dS_dZ = deltaS/deltaZ
    dP_dZ = deltaP / deltaZ
    dMom_dz = deltaMom / deltaZ
    TdS_dZ = T[i]*dS_dZ
    Tdsdz.append(TdS_dZ)
    dpdz.append(dP_dZ)
    outputFile.write(str(zrange[i])+"\t"+str(rho[i])+'\t'+str(-BF_z[i]/rho[i])+"\t"+str(TdS_dZ)+"\t"+str(dP_dZ)+"\t"+str(dMom_dz)+"\n")
        
outputFile.close()