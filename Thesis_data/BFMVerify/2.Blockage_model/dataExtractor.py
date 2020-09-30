# trace generated using paraview version 5.8.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
import sys
import os


R = 287.15
gamma = 1.4

Cp = gamma * R / (gamma - 1)
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
flow = LegacyVTKReader(FileNames=['sym_stator_BFM.vtk'])

z_inlet = -0.09
z_outlet = 0.14

slice1 = Slice(Input=flow)
slice1.SliceType.Origin=[0, 0, z_inlet]
slice1.SliceType.Normal=[0, 0, 1]

IV = IntegrateVariables(Input=slice1)
iv_data = paraview.servermanager.Fetch(IV)

Area = iv_data.GetCellData().GetArray('Area').GetValue(0)
P_in = iv_data.GetPointData().GetArray('Pressure').GetValue(0)/Area
T_in = iv_data.GetPointData().GetArray('Pressure').GetValue(0)/Area

Delete(IV)
weightedP = Calculator(Input=flow)
weightedP.ResultArrayName='weightedP'
weightedP.Function='Momentum_Z*Pressure'

weightedT = Calculator(Input=flow)
weightedT.ResultArrayName='weightedT'
weightedT.Function='Momentum_Z*Temperature'

weightedPtot = Calculator(Input=flow)
weightedPtot.ResultArrayName='weightedPtot'
weightedPtot.Function='Momentum_Z * (Pressure + 0.5 * Density * ((Momentum_X/Density)^2 + (Momentum_Y/Density)^2 + (Momentum_Z/Density)^2))'

weightedTtot = Calculator(Input=flow)
weightedTtot.ResultArrayName='weightedTtot'
weightedTtot.Function='Momentum_Z * (Temperature+ 0.5 * (1 / '+str(Cp)+') * ((Momentum_X/Density)^2 + (Momentum_Y/Density)^2 + (Momentum_Z/Density)^2))'

weightedEntro = Calculator(Input=flow)
weightedEntro.ResultArrayName='weightedEntro'
weightedEntro.Function='Momentum_Z*('+str(Cp)+'*ln(Temperature/'+str(T_in)+') - '+str(R)+'*ln(Pressure/'+str(P_in)+'))'

weightedAngle = Calculator(Input=flow)
weightedAngle.ResultArrayName='weightedAngle'
weightedAngle.Function='Momentum_Z*atan((Momentum_Y*coordsX/sqrt(coordsX^2 + coordsY^2) - Momentum_X*coordsY/sqrt(coordsX^2 + coordsY^2))/Momentum_Z)'

weightedM = Calculator(Input=flow)
weightedM.ResultArrayName='weightedM'
weightedM.Function = 'Momentum_Z*Mach'

axial_range = np.linspace(z_inlet, z_outlet, 1000)

flow2 = AppendAttributes(Input=[flow, weightedP, weightedT, weightedPtot, weightedTtot, weightedEntro, weightedAngle, weightedM])

slice2 = Slice(Input=flow2)
slice2.SliceType.Origin = [0, 0, z_inlet]
slice2.SliceType.Normal = [0, 0, 1]

IV = IntegrateVariables(Input=slice2)
P = []
outfile = open("axialData_BFM.txt", "w+")
outfile.write("#axial coordinate,\tp,\tT,\tp_t,\tT_t,\tds,\talpha,\tMach,\tflux\n")
for z in axial_range:
    print(z)
    slice2.SliceType.Origin = [0, 0, z]
    iv_data = paraview.servermanager.Fetch(IV)
    area = iv_data.GetCellData().GetArray('Area').GetValue(0)

    mom = iv_data.GetPointData().GetArray('Momentum').GetValue(2)
    p = iv_data.GetPointData().GetArray('weightedP').GetValue(0)/mom
    P.append(P)
    T = iv_data.GetPointData().GetArray('weightedT').GetValue(0)/mom
    p_t = iv_data.GetPointData().GetArray('weightedPtot').GetValue(0)/mom
    T_t = iv_data.GetPointData().GetArray('weightedTtot').GetValue(0)/mom
    ds = iv_data.GetPointData().GetArray('weightedEntro').GetValue(0)/mom
    alpha = iv_data.GetPointData().GetArray('weightedAngle').GetValue(0)/mom
    mach = iv_data.GetPointData().GetArray('weightedM').GetValue(0)/mom

    outfile.write("%+16e\t%+16e\t%+16e\t%+16e\t%+16e\t%+16e\t%+16e\t%+16e\t%+16e\n" % (z, p, T, p_t, T_t, ds, alpha, mach, mom/area))

outfile.close()

