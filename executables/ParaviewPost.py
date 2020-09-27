#!/home/evert/anaconda3/bin/python3

import numpy as np
from paraview.simple import *
import os
import sys

HOME = os.environ["M2BFM"]
sys.path.append(HOME + "executables/")

inFile = sys.argv[-1]
from SU2Writer import ReadUserInput
from StagenReader import StagenReader

IN = ReadUserInput(inFile)

class extractData_BFM:
    fileName = "flow_BFM.vtk"
    axial_bounds = []
    radial_bounds = []
    output_datasets = []
    output_data = type('', (), {})()
    fluid_parameters = type('', (), {})()
    calculators = type('', (), {})()
    h = 1e-5

    def __init__(self, IN):
        self.IN = IN
        self.fluid_parameters.gamma = self.IN["gamma"][0]
        self.fluid_parameters.R_gas = self.IN["R_gas"][0]
        self.fluid_parameters.C_p = self.IN["gamma"][0]*self.IN["R_gas"][0]/(self.IN["gamma"][0] - 1)

        os.system("mkdir Performance_Data")
        self.flow_3D_BFM = LegacyVTKReader(FileNames=[self.fileName])
        self.get_row_geom()
        self.get_domain_bounds()
        self.set_functions()
        self.collect_data()
        self.compute_objectives()
    def collect_data(self):
        self.output_datasets = ['X', 'mdot'] + list(self.calculators.__dict__.keys())
        self.extract_axial_data()
        self.write_axial_data()

        self.extract_radial_data()
        self.write_radial_data()
    def compute_objectives(self):
        mdot = self.output_data.axial_data.mdot[0]
        Tt_in = self.output_data.axial_data.wTt[0]
        Tt_out = self.output_data.axial_data.wTt[-1]
        Pt_in = self.output_data.axial_data.wPt[0]
        Pt_out = self.output_data.axial_data.wPt[-1]

        gamma = self.fluid_parameters.gamma
        Tt_out_s = Tt_in * (Pt_out/Pt_in) ** ((gamma - 1)/gamma)

        if self.IN["TYPE"] == 'T':
            eta_tt = (Tt_out - Tt_in)/(Tt_out_s - Tt_in)
        else:
            eta_tt = (Tt_out_s - Tt_in) / (Tt_out - Tt_in)

        power_in = self.fluid_parameters.C_p * mdot * (Tt_out - Tt_in)
        angle_out = self.output_data.axial_data.Alpha[-1]

        self.objective_file = open("Performance_Data/Machine_objectives.txt", "w+")
        self.objective_file.write("total-to-total_efficiency = "+str(eta_tt)+"\n")
        self.objective_file.write("power = "+str(power_in)+"\n")
        self.objective_file.write("mass_flow_rate = " +str(mdot)+"\n")
        self.objective_file.write("outlet_flow_angle = " + str(angle_out) + "\n")
        self.objective_file.close()

    def write_axial_data(self):
        self.output_file = open("Performance_Data/axial_data.txt", "w+")
        self.output_file.write("\t".join([s for s in self.output_datasets])+"\n")
        axial_range = self.output_data.axial_data.X
        for i in range(len(axial_range)):
            section_data = []
            for j in range(len(self.output_datasets)):
                section_data.append(getattr(self.output_data.axial_data, self.output_datasets[j])[i])
            self.output_file.write("\t".join([str(j) for j in section_data])+"\n")
        self.output_file.close()
    def write_radial_data(self):
        station_names = list(self.output_data.radial_data.__dict__.keys())
        for name in station_names:
            rad_file = open("Performance_Data/"+name +".txt", "w+")
            variable_names = list(getattr(self.output_data.radial_data, name).__dict__.keys())
            rad_object = getattr(self.output_data.radial_data, name)
            rad_file.write("\t".join(variable_names)+"\n")
            n = len(rad_object.X)
            for i in range(n):
                rad_data = []
                for var_name in variable_names:
                    rad_data.append(getattr(rad_object, var_name)[i])

                rad_file.write("\t".join([str(s) for s in rad_data])+"\n")
            rad_file.close()

    def extract_radial_data(self):
        radial_range = np.linspace(self.radial_bounds[0], self.radial_bounds[-1], 300)
        stations = [self.axial_bounds[0] + self.h]
        station_names = ['Inlet']
        for i in range(np.shape(self.X_LE)[1]-1):
            stations.append(0.5*(0.5*(self.X_TE[0, i] + self.X_LE[0, i+1]) + 0.5*(self.X_TE[-1, i] + self.X_LE[-1, i+1])))
            station_names.append('Rowgap_'+str(i+1))
        stations.append(self.axial_bounds[-1] - self.h)
        station_names.append('Outlet')
        # stations = [self.axial_bounds[0] + self.h, self.axial_bounds[-1] - self.h]
        # station_names = ['Inlet', 'Outlet']

        slice_ax = Slice(Input=self.flow_data)
        slice_ax.SliceType.Normal = [0, 0, 1]
        slice_ax.SliceType.Origin = [0, 0, stations[0]]

        slice_rad = Slice(Input=slice_ax)
        slice_rad.SliceType = 'Cylinder'
        slice_rad.SliceType.Center = [0.0, 0.0, self.axial_bounds[0]]
        slice_rad.SliceType.Axis = [0.0, 0.0, 1.0]
        slice_rad.SliceType.Radius = radial_range[0]

        IV = IntegrateVariables(Input=slice_rad)
        setattr(self.output_data, 'radial_data', type('', (), {})())
        for i in range(len(stations)):
            setattr(self.output_data.radial_data, station_names[i], type('', (), {})())
            station_object = getattr(self.output_data.radial_data, station_names[i])
            for s in self.output_datasets:
                setattr(station_object, s, [])
            slice_ax.SliceType.Origin = [0, 0, stations[i]]

            for j in range(len(radial_range)):
                slice_rad.SliceType.Radius = radial_range[j]
                iv_data = paraview.servermanager.Fetch(IV)

                if iv_data.GetCellData().GetArray("Length") == None:
                    L = 1.0
                    mom = float('Nan')
                else:
                    L = iv_data.GetCellData().GetArray("Length").GetValue(0)
                    mom = iv_data.GetPointData().GetArray("Momentum").GetValue(2)
                station_object.X.append(radial_range[j])
                station_object.mdot.append(mom/L)

                for s in self.output_datasets[2:]:
                    flow_dat = iv_data.GetPointData().GetArray(s).GetValue(0) / mom
                    output_data_list = getattr(station_object, s)
                    output_data_list.append(flow_dat)
    def extract_axial_data(self):
        axial_range = np.linspace(self.axial_bounds[0] + self.h, self.axial_bounds[1] - self.h, 300)
        setattr(self.output_data, 'axial_data', type('', (), {})())

        for s in self.output_datasets:
            setattr(self.output_data.axial_data, s, [])

        rad_slice = Slice(Input=self.flow_data)
        rad_slice.SliceType.Normal = [0, 0, 1]
        rad_slice.SliceType.Origin = [0, 0, self.axial_bounds[0]+self.h]

        IV = IntegrateVariables(Input=rad_slice)

        for z in axial_range:
            rad_slice.SliceType.Origin = [0, 0, z]
            iv_data = paraview.servermanager.Fetch(IV)
            mom = iv_data.GetPointData().GetArray("Momentum").GetValue(2)
            area = iv_data.GetCellData().GetArray("Area").GetValue(0)

            self.output_data.axial_data.X.append(z)
            self.output_data.axial_data.mdot.append(mom*360/self.IN["WEDGE"][0])

            for s in self.output_datasets[2:]:
                flow_dat = iv_data.GetPointData().GetArray(s).GetValue(0)/mom
                output_data_list = getattr(self.output_data.axial_data, s)
                output_data_list.append(flow_dat)

    def set_functions(self):

        self.define_function(name="wP", input_data=self.flow_3D_BFM, function='Momentum_Z*Pressure')
        self.define_function(name="wT", input_data=self.flow_3D_BFM, function='Momentum_Z*Temperature')

        self.flow_data = AppendAttributes(Input=[getattr(self.calculators, s) for s in list(self.calculators.__dict__.keys())])
        slice_inlet = Slice(Input=self.flow_data)
        slice_inlet.SliceType.Origin = [0, 0, self.axial_bounds[0] + self.h]
        slice_inlet.SliceType.Normal = [0, 0, 1]

        IV = IntegrateVariables(Input=slice_inlet)
        iv_data = paraview.servermanager.Fetch(IV)
        mom_inlet = iv_data.GetPointData().GetArray("Momentum").GetValue(2)
        T_in = iv_data.GetPointData().GetArray("wT").GetValue(0)/mom_inlet
        P_in = iv_data.GetPointData().GetArray("wP").GetValue(0) / mom_inlet

        Delete(IV)
        Delete(slice_inlet)
        Delete(self.flow_data)

        self.define_function(name='ds', input_data=self.flow_3D_BFM,
                             function='Momentum_Z*('+str(self.fluid_parameters.C_p)+'*ln(Temperature/'+str(T_in)+') - '
                                      + str(self.fluid_parameters.R_gas)+'*ln(Pressure/'+str(P_in)+'))')

        self.define_function(name='wPt', input_data=self.flow_3D_BFM,
                             function='Momentum_Z*(Pressure + 0.5 * Density * ((Momentum_X/Density)^2 + '
                                      '(Momentum_Y/Density)^2 + (Momentum_Z/Density)^2))')

        self.define_function(name='wTt', input_data=self.flow_3D_BFM,
                             function='Momentum_Z*(Temperature + 0.5 * (1/'
                                      + str(self.fluid_parameters.C_p)+') * ((Momentum_X/Density)^2 + (Momentum_Y/Density)^2 +(Momentum_Z/Density)^2))')

        self.define_function(name='Alpha', input_data=self.flow_3D_BFM,
                             function='Momentum_Z*(180/'+str(np.pi)+')*atan((Momentum_Y * coordsX/sqrt(coordsX^2 + coordsY^2) - Momentum_X*coordsY/sqrt(coordsX^2 + coordsY^2))/Momentum_Z)')

        self.flow_data = AppendAttributes(Input=[getattr(self.calculators, s) for s in list(self.calculators.__dict__.keys())])

    def define_function(self, name, input_data, function):
        C = Calculator(Input=input_data)
        C.ResultArrayName = name
        C.Function = function
        RenameSource(name, C)
        setattr(self.calculators, name, C)

    def get_domain_bounds(self):

        coords = Calculator(Input=self.flow_3D_BFM)
        coords.ResultArrayName = 'coordinates'
        coords.Function = 'coords'

        boundsfunction = ProgrammableFilter(Input=coords)
        boundsfunction.OutputDataSetType = 'vtkTable'
        boundsfunction.Script = """import numpy as np 
coords = inputs[0].PointData["coordinates"]
    
x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

x_min, y_min, z_min = min(x), min(y), min(z)
x_max, y_max, z_max = max(x), max(y), max(z)

x_bounds = np.array([x_min, x_max])
y_bounds = np.array([y_min, y_max])
z_bounds = np.array([z_min, z_max])
output.RowData.append(x_bounds, "x_bounds")
output.RowData.append(y_bounds, "y_bounds")
output.RowData.append(z_bounds, "z_bounds")"""

        bounds_data = paraview.servermanager.Fetch(boundsfunction)
        self.axial_bounds = [bounds_data.GetRowData().GetArray("z_bounds").GetValue(0),
                        bounds_data.GetRowData().GetArray("z_bounds").GetValue(1)]
        x_bounds = [bounds_data.GetRowData().GetArray("x_bounds").GetValue(0),
                        bounds_data.GetRowData().GetArray("x_bounds").GetValue(1)]
        y_bounds = [bounds_data.GetRowData().GetArray("y_bounds").GetValue(0),
                        bounds_data.GetRowData().GetArray("y_bounds").GetValue(1)]
        self.radial_bounds = [np.sqrt(x_bounds[0]**2 + y_bounds[0]**2), np.sqrt(x_bounds[1]**2 + y_bounds[1]**2)]
        Delete(boundsfunction)
        Delete(coords)

    def get_row_geom(self):
        dir = os.getcwd()
        os.chdir(dir+"/MeangenOutput/")
        S = StagenReader()
        self.X_LE = S.X_LE
        self.X_TE = S.X_TE
        self.Z_LE = S.Z_LE
        self.Z_TE = S.Z_TE
        os.chdir(dir)

extractData_BFM(IN)


# weightedP = Calculator(Input=flow_3D_BFMvtk)
# weightedP.ResultArrayName = 'weightedP'
# weightedP.Function = 'Momentum_Z * Pressure'
#
# weightedT = Calculator(Input=flow_3D_BFMvtk)
# weightedT.ResultArrayName = 'weightedT'
# weightedT.Function = 'Momentum_Z * Temperature'
#
# weightedPtot = Calculator(Input=flow_3D_BFMvtk)
# weightedPtot.ResultArrayName = 'weightedPtot'
# weightedPtot.Function = 'Momentum_Z * (Pressure + 0.5 * Density * ((Momentum_X/Density)^2 + (Momentum_Y/Density)^2 + (Momentum_Z/Density)^2))'
#
# weightedTtot = Calculator(Input=flow_3D_BFMvtk)
# weightedTtot.ResultArrayName = 'weightedTtot'
# weightedTtot.Function = 'Momentum_Z * (Temperature + 0.5 * (1/'+str(Cp)+') * ((Momentum_X/Density)^2 + (Momentum_Y/Density)^2 + (Momentum_Z/Density)^2))'
#
# weightedAngle = Calculator(Input=flow_3D_BFMvtk)
# weightedAngle.ResultArrayName = 'weightedAngle'
# weightedAngle.Function = '-(Momentum_Z)*(180/3.145)*atan((Momentum_Y*(coordsX/sqrt(coordsY^2 + coordsX^2)) - Momentum_X*(coordsY/sqrt(coordsY^2 + coordsX^2)))/Momentum_Z)'
#
# p = ProgrammableFilter(Input=coords)
# p.Script= 'import numpy as np\n' \
#           'C = inputs[0].PointData["coordinates"]\n' \
#           'X = C[:, 0]\n' \
#           'Y = C[:, 1]\n' \
#           'Z = C[:, 2]\n' \
#           'minZ = min(Z)*np.ones(len(X))\n' \
#           'maxZ = max(Z)*np.ones(len(X))\n' \
#           'output.PointData.append(minZ, "minZ")\n' \
#           'output.PointData.append(maxZ, "maxZ")'
#
# P = paraview.servermanager.Fetch(p)
# min_Z = P.GetPointData().GetArray("minZ").GetValue(0)
# max_Z = P.GetPointData().GetArray("maxZ").GetValue(0)
# h = 1e-5
#
# axial_range = np.linspace(min_Z + h, max_Z - h, 200)
#
# flows = AppendAttributes(Input=[weightedP, weightedT])
#
# inlet_slice = Slice(Input=flows)
# inlet_slice.SliceType = 'Plane'
# inlet_slice.SliceType.Origin = [0, 0, axial_range[0]]
# inlet_slice.SliceType.Normal = [0, 0, 1]
#
# integrateVariables1 = IntegrateVariables(Input=inlet_slice)
# iv_data = paraview.servermanager.Fetch(integrateVariables1)
# mom_z = iv_data.GetPointData().GetArray("Momentum").GetValue(2)
# P_s_inlet = iv_data.GetPointData().GetArray("weightedP").GetValue(0)/mom_z
# T_s_inlet = iv_data.GetPointData().GetArray("weightedT").GetValue(0)/mom_z
#
# weightedEntro = Calculator(Input = flow_3D_BFMvtk)
# weightedEntro.ResultArrayName = 'weightedEntro'
# weightedEntro.Function = 'Momentum_Z * ('+str(Cp)+'*ln(Temperature / '+str(T_s_inlet) + ') - 287.14*ln(Pressure / '+str(P_s_inlet)+'))'
#
# massrate = Calculator(Input=flow_3D_BFMvtk)
# massrate.ResultArrayName = 'massrate'
# massrate.Function = 'Momentum_Z * Blockage_Factor'
#
# Delete(flows)
# Delete(inlet_slice)
# flowvars = AppendAttributes(Input=[weightedP, weightedT, weightedPtot, weightedTtot, weightedAngle, weightedEntro, massrate])
#
# slice1 = Slice(Input = flowvars)
# slice1.SliceType = 'Plane'
# slice1.SliceType.Origin = [0, 0, axial_range[0]]
# slice1.SliceType.Normal = [0, 0, 1]
#
#
#
# outfile = open("AxialData.txt", "w+")
# enthalpy = []
# T = []
# P = []
# Ttot = []
# Ptot = []
# mflux = []
# angle = []
# for i in range(len(axial_range)):
#     slice1.SliceType.Origin = [0, 0, axial_range[i]]
#     IV = IntegrateVariables(Input=slice1)
#     iv_data = paraview.servermanager.Fetch(IV)
#     momentum = iv_data.GetPointData().GetArray("Momentum").GetValue(2)
#     area = iv_data.GetCellData().GetArray("Area").GetValue(0)
#     T_stat = iv_data.GetPointData().GetArray("weightedT").GetValue(0)/momentum
#     P_stat = iv_data.GetPointData().GetArray("weightedP").GetValue(0)/momentum
#     T_tot = iv_data.GetPointData().GetArray("weightedTtot").GetValue(0)/momentum
#     P_tot = iv_data.GetPointData().GetArray("weightedPtot").GetValue(0) / momentum
#     alpha = iv_data.GetPointData().GetArray("weightedAngle").GetValue(0) / momentum
#     mrate = iv_data.GetPointData().GetArray("massrate").GetValue(0) / area
#     enthalpy.append(momentum*360*Cp*T_tot)
#     T.append(T_stat)
#     P.append(P_stat)
#     Ttot.append(T_tot)
#     Ptot.append(P_tot)
#     mflux.append(momentum)
#     angle.append(alpha)
#     outfile.write("%+16e\t%+16e\t%+16e\t%+16e\t%+16e\n" % (axial_range[i], T_stat, P_stat, T_tot, P_tot))
#
# outfile.close()
#
# w_isentropic = Cp * Ttot[0]*(1 - (Ptot[-1]/Ptot[0])**((gamma - 1)/gamma))
# w_real = Cp * (Ttot[0] - Ttot[-1])
# eta = w_real / w_isentropic
# power = Cp*360*(mflux[-1]*Ttot[-1] - mflux[0]*Ttot[0])
# print(w_isentropic)
# print(w_real)
# print("Outlet flow angle: " + str(angle[-1]) + " deg")
