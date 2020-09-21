import numpy as np
from paraview.simple import *
import os
Cp = 1150
gamma = 1.33

fileName = "flow_BFM.vtk"
flow_3D_BFMvtk = LegacyVTKReader(FileNames=[fileName])

coords = Calculator(Input=flow_3D_BFMvtk)
coords.ResultArrayName = 'coordinates'
coords.Function = 'coords'

weightedP = Calculator(Input=flow_3D_BFMvtk)
weightedP.ResultArrayName = 'weightedP'
weightedP.Function = 'Momentum_Z * Pressure'

weightedT = Calculator(Input=flow_3D_BFMvtk)
weightedT.ResultArrayName = 'weightedT'
weightedT.Function = 'Momentum_Z * Temperature'

weightedPtot = Calculator(Input=flow_3D_BFMvtk)
weightedPtot.ResultArrayName = 'weightedPtot'
weightedPtot.Function = 'Momentum_Z * (Pressure + 0.5 * Density * ((Momentum_X/Density)^2 + (Momentum_Y/Density)^2 + (Momentum_Z/Density)^2))'

weightedTtot = Calculator(Input=flow_3D_BFMvtk)
weightedTtot.ResultArrayName = 'weightedTtot'
weightedTtot.Function = 'Momentum_Z * (Temperature + 0.5 * (1/'+str(Cp)+') * ((Momentum_X/Density)^2 + (Momentum_Y/Density)^2 + (Momentum_Z/Density)^2))'

weightedAngle = Calculator(Input=flow_3D_BFMvtk)
weightedAngle.ResultArrayName = 'weightedAngle'
weightedAngle.Function = '-(Momentum_Z)*(180/3.145)*atan((Momentum_Y*(coordsX/sqrt(coordsY^2 + coordsX^2)) - Momentum_X*(coordsY/sqrt(coordsY^2 + coordsX^2)))/Momentum_Z)'

p = ProgrammableFilter(Input=coords)
p.Script= 'import numpy as np\n' \
          'C = inputs[0].PointData["coordinates"]\n' \
          'X = C[:, 0]\n' \
          'Y = C[:, 1]\n' \
          'Z = C[:, 2]\n' \
          'minZ = min(Z)*np.ones(len(X))\n' \
          'maxZ = max(Z)*np.ones(len(X))\n' \
          'output.PointData.append(minZ, "minZ")\n' \
          'output.PointData.append(maxZ, "maxZ")'

P = paraview.servermanager.Fetch(p)
min_Z = P.GetPointData().GetArray("minZ").GetValue(0)
max_Z = P.GetPointData().GetArray("maxZ").GetValue(0)
h = 1e-5

axial_range = np.linspace(min_Z + h, max_Z - h, 200)

flows = AppendAttributes(Input=[weightedP, weightedT])

inlet_slice = Slice(Input=flows)
inlet_slice.SliceType = 'Plane'
inlet_slice.SliceType.Origin = [0, 0, axial_range[0]]
inlet_slice.SliceType.Normal = [0, 0, 1]

integrateVariables1 = IntegrateVariables(Input=inlet_slice)
iv_data = paraview.servermanager.Fetch(integrateVariables1)
mom_z = iv_data.GetPointData().GetArray("Momentum").GetValue(2)
P_s_inlet = iv_data.GetPointData().GetArray("weightedP").GetValue(0)/mom_z
T_s_inlet = iv_data.GetPointData().GetArray("weightedT").GetValue(0)/mom_z

weightedEntro = Calculator(Input = flow_3D_BFMvtk)
weightedEntro.ResultArrayName = 'weightedEntro'
weightedEntro.Function = 'Momentum_Z * ('+str(Cp)+'*ln(Temperature / '+str(T_s_inlet) + ') - 287.14*ln(Pressure / '+str(P_s_inlet)+'))'

massrate = Calculator(Input=flow_3D_BFMvtk)
massrate.ResultArrayName = 'massrate'
massrate.Function = 'Momentum_Z * Blockage_Factor'

Delete(flows)
Delete(inlet_slice)
flowvars = AppendAttributes(Input=[weightedP, weightedT, weightedPtot, weightedTtot, weightedAngle, weightedEntro, massrate])

slice1 = Slice(Input = flowvars)
slice1.SliceType = 'Plane'
slice1.SliceType.Origin = [0, 0, axial_range[0]]
slice1.SliceType.Normal = [0, 0, 1]



outfile = open("AxialData.txt", "w+")
enthalpy = []
T = []
P = []
Ttot = []
Ptot = []
mflux = []
angle = []
for i in range(len(axial_range)):
    slice1.SliceType.Origin = [0, 0, axial_range[i]]
    IV = IntegrateVariables(Input=slice1)
    iv_data = paraview.servermanager.Fetch(IV)
    momentum = iv_data.GetPointData().GetArray("Momentum").GetValue(2)
    area = iv_data.GetCellData().GetArray("Area").GetValue(0)
    T_stat = iv_data.GetPointData().GetArray("weightedT").GetValue(0)/momentum
    P_stat = iv_data.GetPointData().GetArray("weightedP").GetValue(0)/momentum
    T_tot = iv_data.GetPointData().GetArray("weightedTtot").GetValue(0)/momentum
    P_tot = iv_data.GetPointData().GetArray("weightedPtot").GetValue(0) / momentum
    alpha = iv_data.GetPointData().GetArray("weightedAngle").GetValue(0) / momentum
    mrate = iv_data.GetPointData().GetArray("massrate").GetValue(0) / area
    enthalpy.append(momentum*360*Cp*T_tot)
    T.append(T_stat)
    P.append(P_stat)
    Ttot.append(T_tot)
    Ptot.append(P_tot)
    mflux.append(momentum)
    angle.append(alpha)
    outfile.write("%+16e\t%+16e\t%+16e\t%+16e\t%+16e\n" % (axial_range[i], T_stat, P_stat, T_tot, P_tot))

outfile.close()

w_isentropic = Cp * Ttot[0]*(1 - (Ptot[-1]/Ptot[0])**((gamma - 1)/gamma))
w_real = Cp * (Ttot[0] - Ttot[-1])
eta = w_real / w_isentropic
power = Cp*360*(mflux[-1]*Ttot[-1] - mflux[0]*Ttot[0])
print(w_isentropic)
print(w_real)
print("Outlet flow angle: " + str(angle[-1]) + " deg")
