# --------------------------------------------------------------------------------------------------------------- #
# This set of classes is used for creating a 2D mesh for BFM or physical blade analysis purposes. It does so through the
# use ofUMG2. In the future, this might be changed to Gmesh, due to the latter being somewhat faster and more easy to
# use.
# For the class to work, UMG2 has to have been downloaded and its installation location has to be appended to the
# PATH according to installation instructions.
# The WriteUMG class is used for writing UMG2 files for each individual mesh zone for either the BFM mesh or physical
# blade mesh. The writeStageMesh class then combines these meshes into a single mesh for the entire machine.
# --------------------------------------------------------------------------------------------------------------- #
import numpy as np
from scipy.interpolate import interp1d
import os

class WriteUMG:
    def __init__(self, rowNumber, stage, Meangen, IN, bodyForce, blade):

        self.n_blade = 1            # Number blades in the mesh.
        self.rowNumber = rowNumber  # Row number.
        self.stage = stage          # Stage number.
        self.Meangen = Meangen      # Storing Meangen class.
        self.name = 'turbo'         # Mesh name.

        thisDir = os.getcwd()       # Getting current directory.

        self.wedge = IN["WEDGE"][0]

        # Getting the number of axial points for the blade passage, boundary layer count and boundary layer thickness.
        np = IN["AXIAL_POINTS"][0]
        self.n_bl = int(IN["BOUNDARY_LAYER_COUNT"][0])
        self.thc_bl = IN["BOUNDARY_LAYER_THICKNESS"][0]

        # Determining the axial chord length. In case of a compressor, the rotor chord is followed by the stator chord.
        # In case of a turbine, the stator is followed by the rotor.
        if self.Meangen.machineType == 'C':
            if rowNumber == 1:
                chord = self.Meangen.chord_R[stage-1]
            else:
                chord = self.Meangen.chord_S[stage-1]
        else:
            if rowNumber == 1:
                chord = self.Meangen.chord_S[stage-1]
            else:
                chord = self.Meangen.chord_R[stage-1]

        # Calculating the inflow cell sizes. These are the passage cell sizes multiplied by the specified inlet factor.
        if self.rowNumber == 1 and self.stage == 1:
            self.h_max_inflow = IN["INLET_FACTOR"][0] * chord/np
            self.h_min_inflow = self.h_max_inflow / IN["INLET_FACTOR"][0]
        else:
            self.h_max_inflow = chord / np
            self.h_min_inflow = self.h_max_inflow
        self.radCrv_inflow = 5

        # Calculating the passage cell sizes.
        self.h_max_perio = chord/np
        self.h_min_perio = self.h_max_perio
        self.radCrv_perio = 5

        # Calculating the blade surface cell sizes.
        self.h_max_blade = IN["BLADE_FACTOR"][0] * self.h_max_perio
        self.h_min_blade = self.h_max_blade / 10
        self.radCrv_blade = 5

        # Calculating the outflow cell sizes. These are the passage cell sizes multiplied by the specified outlet factor
        if self.rowNumber == 2:
            self.h_max_outflow = IN["OUTLET_FACTOR"][0] * chord/np
            self.h_min_outflow = self.h_max_outflow / IN["OUTLET_FACTOR"][0]
        else:
            self.h_max_outflow = chord / np
            self.h_min_outflow = self.h_max_outflow
        self.radCrv_outflow = 5

        # Getting installation directory.
        HOME = os.environ["M2BFM"]

        # Copying createmesh template file from templates directory.
        os.system("cp " + HOME + "templates/createmesh.template ./")

        # Mesh zone is created according to the user's mesh specification.
        if bodyForce:
            print("Starting BFM mesh computation")
            self.makeBFMMesh()
        else:
            print("No BFM mesh creation requested")
        os.chdir(thisDir)
        if blade:
            print("Starting blade mesh computation")
            self.makeBladeMesh()
        else:
            print("No blade mesh creation requested")

    def makeGeom_BFM(self):
        # This function is used for creating the geometry configuration file for the BFM mesh.

        # Getting the leading and trailing edge coordinates from the Meangen class.
        x_le = self.Meangen.X_LE
        x_te = self.Meangen.X_TE

        # Calculating the blade pitch according to the machine type.
        if self.Meangen.machineType == 'C':
            if self.rowNumber % 2 != 0:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1])/self.Meangen.N_b_R[self.stage - 1]
            else:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1]) / self.Meangen.N_b_S[self.stage - 1]
        else:
            if self.rowNumber % 2 != 0:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1])/self.Meangen.N_b_S[self.stage - 1]
            else:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1]) / self.Meangen.N_b_R[self.stage - 1]
        self.Pitch = pitch

        # Setting relevant x coordinates for the mesh geometry.
        x1 = x_le[1, 2 * (self.stage - 1) + self.rowNumber - 1]

        x2 = x1
        y2 = pitch * self.n_blade

        x3 = x_te[1, 2 * (self.stage - 1) + self.rowNumber - 1]

        x4 = x3

        # Depending on the row and stage number, either a far field, outflow or row or stage gap mesh zone is created.
        if self.stage == 1 and self.rowNumber == 1:
            x_fwd = x1 - (x3 - x1)
            x_bck = 0.5 * (x_te[1, self.rowNumber - 1] + self.Meangen.X_LE[1, self.rowNumber])
        elif self.stage == self.Meangen.n_stage and self.rowNumber == 2:
            x_fwd = 0.5 * (x_le[1, 2 * (self.stage - 1) + self.rowNumber - 1] + self.Meangen.X_TE[1, 2 * (self.stage - 1) + self.rowNumber - 2])
            x_bck = x3 + (x3 - x1)
        else:
            x_fwd = 0.5 * (x_le[1, 2 * (self.stage - 1) + self.rowNumber - 1] + self.Meangen.X_TE[1, 2 * (self.stage - 1) + self.rowNumber - 2])
            x_bck = 0.5 * (x_te[1, 2 * (self.stage - 1) + self.rowNumber - 1] + self.Meangen.X_LE[1, 2 * (self.stage - 1) + self.rowNumber])

        # Collecting all x-coordinates of the mesh corners.
        self.X_curve = [[x_fwd, x1], [x1, x1], [x1, x_fwd], [x_fwd, x_fwd],
                        [x1, x4], [x4, x4], [x3, x2], [x2, x1],
                        [x4, x_bck], [x_bck, x_bck], [x_bck, x3], [x3, x4]]

        # Collecting all y-coordinates of the mesh corners.
        self.Y_curve = [[0, 0], [0, y2], [y2, y2], [y2, 0],
                        [0, 0], [0, y2], [y2, y2], [y2, 0],
                        [0, 0], [0, y2], [y2, y2], [y2, 0]]

        # Setting names to the mesh boundaries.
        self.names = ["PERIO_DOWN", "OUTFLOW", "PERIO_UP", "INFLOW",
                      "PERIO_DOWN", "OUTFLOW", "PERIO_UP", "INFLOW",
                      "PERIO_DOWN", "OUTFLOW", "PERIO_UP", "INFLOW"]

        # Specifying boundary types.
        self.types = [4, 3, 5, 1,
                      4, 3, 5, 1,
                      4, 3, 5, 1]

        # Specifying boundary periodicity.
        self.periodic = [3, 0, 1, 0,
                         3, 0, 1, 0,
                         3, 0, 1, 0]

        # Specifing the boundary order.
        self.order = [1, 2, 3, 4,
                      1, 2, 3, 4,
                      1, 2, 3, 4]

    def makeGeom_blade(self):
        # This function is used for generating the mesh geometry configuration file for physical blade analysis.

        # Getting the leading and trailing edge coordinates from the Meangen class.
        x_le = self.Meangen.X_LE
        x_te = self.Meangen.X_TE

        # Getting the blade surface coordinates from Parablade output.
        coordDir = os.getcwd() + "/output/coordinates/surface_coordinates.csv"
        p, x, y, u, v = np.loadtxt(coordDir, unpack=True, delimiter=',\t', skiprows=1)

        # Constructing upper and lower blade surface curves.
        i_min = list(x).index(min(x))
        i_max = list(x).index(max(x))
        if i_min < i_max:
            x_down = x[i_min:i_max+1]
            y_down = y[i_min:i_max+1]
            x_up = np.concatenate((x[i_max:], x[:i_min+1]), axis=0)
            y_up = np.concatenate((y[i_max:], y[:i_min + 1]), axis=0)
        else:
            x_down = np.concatenate((x[i_min:], x[:i_max+1]), axis=0)
            y_down = np.concatenate((y[i_min:], y[:i_max+1]), axis=0)
            x_up = x[i_max:i_min+1]
            y_up = y[i_max:i_min+1]

        # Setting up interpolation function for upper blade surface curve.
        Y_up = interp1d(x_up, y_up, kind='linear')

        # Interpolating upper blade surface curve at lower blade surface curve x-coordinates.
        y_up = Y_up(x_down)

        # Setting up a curve which is the average of the upper and lower blade surface curve. This will be used for the
        # periodic boundaries.
        x_av = x_down
        y_av = 0.5*(y_up + y_down)

        # Determining leading and trailing edge coordinates.
        x_0 = x_av[0]
        y_0 = y_av[0]
        x_1 = x_av[-1]
        y_1 = y_av[-1]

        # Calculating blade pitch according to machine type.
        if self.Meangen.machineType == 'C':
            if self.rowNumber % 2 != 0:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1])/self.Meangen.N_b_R[self.stage - 1]
            else:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1]) / self.Meangen.N_b_S[self.stage - 1]
        else:
            if self.rowNumber % 2 != 0:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1])/self.Meangen.N_b_S[self.stage - 1]
            else:
                pitch = (2 * np.pi * self.Meangen.r_m[self.stage - 1]) / self.Meangen.N_b_R[self.stage - 1]
        self.Pitch = pitch

        # Depending on row and stage number, a far field, row gap or stage gap is added to the front or rear of the
        # blade passage mesh.
        if self.stage == 1 and self.rowNumber == 1:
            x_2 = x_0 - 2*(x_1 - x_0)
            x_6 = 0.5 * (x_te[1, 2 * (self.stage - 1) + self.rowNumber - 1] + self.Meangen.X_LE[1, 2 * (self.stage - 1) + self.rowNumber])
        elif self.stage == self.Meangen.n_stage and self.rowNumber == 2:
            x_2 = 0.5 * (x_le[1, 2 * (self.stage - 1) + self.rowNumber - 1] + self.Meangen.X_TE[1, 2 * (self.stage - 1) + self.rowNumber - 2])
            x_6 = x_1 + 2*(x_1 - x_0)
        else:
            x_2 = 0.5 * (x_le[1, 2 * (self.stage - 1) + self.rowNumber - 1] + self.Meangen.X_TE[1, 2 * (self.stage - 1) + self.rowNumber - 2])
            x_6 = 0.5 * (x_te[1, 2 * (self.stage - 1) + self.rowNumber - 1] + self.Meangen.X_LE[1, 2 * (self.stage - 1) + self.rowNumber])

        # Setting up the coordinates for the mesh corners.
        y_2 = y_0 - 0.5 * pitch

        x_3 = x_2
        y_3 = y_0 + (0.5 + self.n_blade - 1) * pitch

        x_4 = x_0
        y_4 = y_3

        x_5 = x_1
        y_5 = y_1 + (0.5 + self.n_blade - 1) * pitch

        y_6 = y_5
        x_7 = x_6
        y_7 = y_1 - 0.5 * pitch

        x_8 = x_1
        y_8 = y_7

        x_9 = x_0
        y_9 = y_0 - 0.5 * pitch

        # Specifing boundary names.

        self.names = ["BLADE", "BLADE", "INFLOW", "PER_INF_DOWN", "PER_CHANNEL_DOWN", "PER_GAP_DOWN",
                      "PER_INF_UP", "PER_CHANNEL_UP", "PER_GAP_UP", "OUTFLOW"]

        # Specifying boundary types.

        self.types = [8, 8, 1, 4, 4, 4, 5, 5, 5, 3]

        # Specifying boundary order.
        self.order = [1, 2, 4, 5, 6, 10, 9, 8, 7, 3]

        # Specifying periodic boundary order.
        self.periodic = [0, 0, 0, 7, 8, 9, 4, 5, 6, 0]

        # Setting up curve coordinates.
        self.X_curve = [x_down, x_down, [x_2, x_3], [x_2, x_9], x_av, [x_8, x_7], [x_3, x_4], x_av, [x_5, x_6], [x_7, x_6]]
        self.Y_curve = [y_down, y_up, [y_2, y_3], [y_2, y_9], y_av - 0.5 * pitch, [y_8, y_7], [y_3, y_4],
                        y_av + (0.5 + self.n_blade - 1) * pitch, [y_5, y_6], [y_7, y_6]]


    def makeBFMMesh(self):
        # This function writes the UMG2 configuration files for the construction of the body-force mesh. Once written,
        # UMG2 uses these files to create a mesh for each mesh zone. These mesh files are then appended into a single
        # mesh file for the blade row.

        # Creating mesh geometry curves and coordinates.
        self.makeGeom_BFM()

        # Creating a folder where body-force mesh files are stored.
        if os.path.exists("BFMMesh"):
            print("Directory already exists, rewriting files")
        else:
            os.system("mkdir BFMMesh")
        meshDir = os.getcwd() + "/BFMMesh"
        os.chdir(meshDir)

        # Creating Db directory for UMG2 to read from.
        if not os.path.exists("Db"):
            os.system("mkdir Db")
        fileDir = os.getcwd() + "/Db/"

        # Specifying names for the three zones that make up each blade row mesh.
        zone_names = ["inlet", "channel", "outlet"]

        # Opening the file used for the blade row mesh.
        meshFile = open("BFM_mesh.su2", "w+")
        meshFile.write("NZONE=  "+str(len(zone_names))+"\n")

        # Looping over the three zones in the blade row mesh to create a UMG2 mesh for each zone.
        for k in range(len(zone_names)):
            os.chdir(fileDir)

            # Opening the geometry file to write all geometric properties specified for the body-force mesh.
            geomFile = open("geometry."+zone_names[k], "w+")
            geomFile.write("Number of surfaces\n%i\n" % 4)
            for i in range(4*k, 4*(k + 1)):
                geomFile.write("%s\n' S '\n" % (self.names[i]))
                geomFile.write("dim\tnp\n2\t%i\nx\ty\n" % (len(self.X_curve[i])))
                for j in range(len(self.X_curve[i])):
                    geomFile.write("%+.5e\t%+.5e\n" % (self.X_curve[i][j], self.Y_curve[i][j]))
            geomFile.close()

            # Opening the topology file to specify all boundary types.
            topoFile = open("topology."+zone_names[k], "w+")
            topoFile.write("curve type\tperiodic curve\tModifiable curve\n")
            for i in range(4*k, 4*(k + 1)):
                topoFile.write("%i\t%i\t%i\n" % (self.types[i], self.periodic[i], 0))
            topoFile.write("Number of ZONE\n")
            topoFile.write("1\nZONE 1\n")
            for i in range(4*(k-1), 4*k - 1):
                topoFile.write(" %i\n" % (self.order[i]))
            topoFile.write("%i\n" % (-(self.order[-1])))
            topoFile.close()

            # Opening the spacing control file to write the mesh cell sizes for each boundary.
            spacingFile = open("spacingcontrol."+zone_names[k], "w+")
            spacingFile.write("thk_bl\tn\tBC\tGEOM\tCV\n")
            spacingFile.write("%+.5e\t%i\taxl\t0\n\n" % (self.thc_bl, 5))
            spacingFile.write("PITCH\txc\tyc\n")
            spacingFile.write("%+.5e\t1.0\t1.0\n\n" % (self.Pitch))
            spacingFile.write("1\tINFLOW\th_min\th_max\tNode per RadCRv\n")
            if k == 0:
                spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_max_inflow, self.h_max_inflow, self.radCrv_inflow))
            else:
                spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_inflow, self.h_min_inflow, self.radCrv_inflow))
            spacingFile.write("8\tBLADE\th_min\th_max\tNode per RadCRv\n")
            spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_blade, self.h_max_blade, self.radCrv_blade))
            spacingFile.write("3\tOUTFLOW\th_min\th_max\tNode per RadCRv\n")
            if k == 2:
                spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_max_outflow, self.h_max_outflow, self.radCrv_outflow))
            else:
                spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_outflow, self.h_min_outflow, self.radCrv_outflow))

            spacingFile.write("4\tPERIO\th_min\th_max\tNode per RadCRv\n")
            spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
            spacingFile.write("5\tPERIO\th_min\th_max\tNode per RadCRv\n")
            spacingFile.write("%+.5e\t%+.5e\t%i\n\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
            spacingFile.write("NZONES\n1\n")
            spacingFile.close()

            # Opening the options file.
            optionsFile = open("options", "w+")
            optionsFile.write("fmt\tname\n")
            optionsFile.write("'grd'\t'%s'\n" % zone_names[k])
            optionsFile.write("optimization\n1\nmax element deformation\n1.0\nlayer of background grid\n3\n")
            optionsFile.write("Periodic geometry\n.true\t1.e-6\nScaling for SU2 file\n1.0\nnumber of boundary layers\n")
            optionsFile.write("0\t%+.5e\n" % self.thc_bl)
            optionsFile.write("Graph for hybrid mesh construction\n.true\nKind of radial basis function(1-11)\n11\n")
            optionsFile.write("Support radius for compact basis functions\n0.05\n")
            optionsFile.close()

            # Writing the UMG2 mesh file for the body-force mesh of the current zone.
            os.chdir(meshDir)
            print("Writing body-force mesh for the " + zone_names[k] + " of blade row " + str(self.rowNumber) + "...")
            os.system("HYMESH.sh > UMG2_"+zone_names[k]+".out")
            print("Done!")

            # Copying SU2_PERIO configuration template file to mesh zone directory.
            os.system("cp ../createmesh.template ./createmesh.cfg")
            # Replacing pitch template value.
            os.system("sed -i 's/PITCH/%+.5e/g' createmesh.cfg" % self.Pitch)

            # Executing SU2_PERIO to create a periodic mesh zone.
            os.system("SU2_PERIO < createmesh.cfg > SU2_PERIO_Output")

            # Renaming mesh zone to reflect the location in the stage.
            os.system("mv ./mesh_out.su2 ./mesh_"+zone_names[k]+".su2")
            os.system("mv ./tec_mesh.dat ./mesh_tec_" + zone_names[k] + ".dat")

            # Replacing mesh boundary names to reflect location in the stage.
            os.system("sed -i 's/inflow/inflow_"+str(k+1)+"/g' mesh_"+zone_names[k]+".su2")
            os.system("sed -i 's/outflow/outflow_" + str(k + 1) + "/g' mesh_" + zone_names[k] + ".su2")
            os.system("sed -i 's/periodic1/periodic1_" + str(k + 1) + "/g' mesh_" + zone_names[k] + ".su2")
            os.system("sed -i 's/periodic2/periodic2_" + str(k + 1) + "/g' mesh_" + zone_names[k] + ".su2")
            os.system("sed -i 's|IZONE=  1|IZONE=  "+str(k + 1) + "|' mesh_" + zone_names[k] + ".su2")

            # Opening zone mesh file for reading.
            current_mesh = open("mesh_" + zone_names[k] + ".su2", "r")
            # Skipping first line in the mesh zone file.
            next(current_mesh)
            # Writing
            lines = current_mesh.readlines()[1:]
            meshFile.writelines(lines)
            # for line in current_mesh:
            #     meshFile.write(line)



    def makeBladeMesh(self):
        # This function creates the 2D physical blade mesh.
        self.makeGeom_blade()             # Writing geometry file.
        # Creating a folder where all relevant blade meshing input files are stored.
        if os.path.exists("BladeMesh"):
            print("Directory already exists, rewriting files")
        else:
            os.system("mkdir BladeMesh")

        meshDir = os.getcwd() + "/BladeMesh"
        os.chdir(meshDir)

        # Making Db folder for UMG2 to read from.
        if not os.path.exists("Db"):
            os.system("mkdir Db")
        fileDir = os.getcwd() + "/Db/"
        os.chdir(fileDir)

        # Opening geometry configuration file and writing all blade mesh curves.
        geomFile = open("geometry." + self.name, "w+")
        geomFile.write("Number of surfaces\n%i\n" % (len(self.names)))
        for i in range(len(self.names)):
            geomFile.write("%s\n' S '\n" % (self.names[i]))
            geomFile.write("dim\tnp\n2\t%i\nx\ty\n" % (len(self.X_curve[i])))
            for j in range(len(self.X_curve[i])):
                geomFile.write("%+.5e\t%+.5e\n" % (self.X_curve[i][j], self.Y_curve[i][j]))
        geomFile.close()

        # Opening topology configuration file and writing curve types and order.
        topoFile = open("topology." + self.name, "w+")
        topoFile.write("curve type\tperiodic curve\tModifiable curve\n")
        for i in range(len(self.types)):
            topoFile.write("%i\t%i\t%i\n" % (self.types[i], max([0, self.periodic[i]]), 0))
        topoFile.write("Number of ZONE\n")
        topoFile.write("1\nZONE 1\n")
        for i in range(len(self.order) - 1):
            topoFile.write(" %i\n" % (self.order[i]))
        topoFile.write("%i\n" % (-(self.order[-1])))
        topoFile.close()

        # Opening spacing control file and writing down all mesh size parameters.
        spacingFile = open("spacingcontrol." + self.name, "w+")
        spacingFile.write("thk_bl\tn\tBC\tGEOM\tCV\n")
        spacingFile.write("%+.5e\t%i\taxl\t0\n\n" % (self.thc_bl, self.n_bl))
        spacingFile.write("PITCH\txc\tyc\n")
        spacingFile.write("%+.5e\t1.0\t1.0\n\n" % (self.Pitch))
        spacingFile.write("1\tINFLOW\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_inflow, self.h_max_inflow, self.radCrv_inflow))
        spacingFile.write("8\tBLADE\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_blade, self.h_max_blade, self.radCrv_blade))
        spacingFile.write("3\tOUTFLOW\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_outflow, self.h_max_outflow, self.radCrv_outflow))
        spacingFile.write("4\tPERIO\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
        spacingFile.write("5\tPERIO\th_min\th_max\tNode per RadCRv\n")
        spacingFile.write("%+.5e\t%+.5e\t%i\n\n" % (self.h_min_perio, self.h_max_perio, self.radCrv_perio))
        spacingFile.write("NZONES\n1\n")
        spacingFile.close()

        # Opening and writing the options configuration file.
        optionsFile = open("options", "w+")
        optionsFile.write("fmt\tname\n")
        optionsFile.write("'grd'\t'%s'\n" % (self.name))
        optionsFile.write("optimization\n1\nmax element deformation\n1.0\nlayer of background grid\n3\n")
        optionsFile.write("Periodic geometry\n.true\t1.e-6\nScaling for SU2 file\n1.0\nnumber of boundary layers\n")
        optionsFile.write("%i\t%+.5e\n" % (self.n_bl, self.thc_bl))
        optionsFile.write("Graph for hybrid mesh construction\n.true\nKind of radial basis function(1-11)\n11\n")
        optionsFile.write("Support radius for compact basis functions\n0.05\n")
        optionsFile.close()

        # Executing HYMESH to construct hybrid mesh for blade row.
        os.chdir(meshDir)
        print("Constructing hybrid mesh for blade row number " + str(self.rowNumber) + "...")
        os.system("HYMESH.sh > Blade_Mesh.out")
        print("Done!")

        # Copying SU2_PERIO configuration template file and replacing the pitch template variable with the actual pitch.
        os.system("cp ../createmesh.template ./createmesh.cfg")
        os.system("sed -i 's/PITCH/%+.5e/g' createmesh.cfg" % self.Pitch)

        # Executing SU2_PERIO
        print("Constructing periodic mesh for blade row number " + str(self.rowNumber) + "...")
        os.system("SU2_PERIO < createmesh.cfg > SU2_PERIO_blade.out")
        print("Done!")

class writeStageMesh_BFM:
    # This class combines the different UMG2 mesh zones into a single mesh which can be used for full-machine analysis
    # In SU2 for body-force analysis.
    def __init__(self, Meangen):
        self.Meangen = Meangen  # Storing Meangen class.
        self.dir = os.getcwd()  # Getting current directory.
        self.n_stage = Meangen.n_stage  # Number of stages.
        self.n_rows = self.n_stage * 2  # Number of blade rows.
        self.n_zone = 2 * self.n_stage * 3    # Total number of mesh zones.

        # Replacing boundary names in the individual mesh zone files.
        self.replaceTerms()

        # Opening stage mesh SU2 file.
        self.meshFile = open(self.dir + "/BFM_mesh_machine.su2", "w+")

        # Writing total number of zones.
        self.meshFile.write("NZONE=  %i\n" % self.n_zone)
        self.writeMeshFile()
        self.meshFile.close()
    def replaceTerms(self):
        # This function replaces the boundary names of the individual zone boundaries with names compatible with the
        # full machine mesh.
        k = 0

        # Looping over the stages and blade rows to replace the boundary names in the inlet, channel and outlet zone
        # of each blade row mesh.
        for i in range(self.n_stage):
            for j in [1, 2]:
                meshDir = self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/BFMMesh/"
                os.chdir(meshDir)
                for q in range(1, 4):
                    os.system("sed -i 's|IZONE=  "+str(q)+"|IZONE= "+str(k*3 + q)+"|' BFM_mesh.su2")
                    os.system("sed -i 's/inflow_"+str(q)+"/inflow_" + str(k*3 + q) + "/g' BFM_mesh.su2")
                    os.system("sed -i 's/outflow_"+str(q)+"/outflow_" + str(k * 3 + q) + "/g' BFM_mesh.su2")
                    os.system("sed -i 's/periodic1_" + str(q) + "/periodic1_" + str(k * 3 + q) + "/g' BFM_mesh.su2")
                    os.system("sed -i 's/periodic2_" + str(q) + "/periodic2_" + str(k * 3 + q) + "/g' BFM_mesh.su2")
                k += 1

    def writeMeshFile(self):
        # This function combines the blade row meshes into a single SU2 mesh file which can be used as input for SU2.
        for i in range(self.n_stage):
            for j in [1, 2]:
                with open(self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/BFMMesh/BFM_mesh.su2", "r") as BFMmesh:
                    # Reading the lines in the blade row mesh after the first line and writing them into the machine
                    # mesh.
                    lines = BFMmesh.readlines()[1:]
                    self.meshFile.writelines(lines)
                BFMmesh.close()


class writeStageMesh_Blade:
    # This class combines the physical blade meshes of the machine blade rows into a single SU2 mesh which can be used
    # for 2D analysis in SU2.
    def __init__(self, Meangen):
        self.Meangen = Meangen  # Storing the Meangen class.
        self.dir = os.getcwd()  # Getting current directory.
        self.n_stage = Meangen.n_stage  # Number of stages.
        self.n_rows = self.n_stage * 2  # Number of blade rows.
        self.n_zone = self.n_rows   # Number of mesh zones.

        # Replacing the boundary names with terms compatible with a full-machine mesh.
        self.replaceTerms()

        # Opening machine mesh file and writing total number of zones.
        self.meshFile = open(self.dir + "/Blade_mesh_machine.su2", "w+")
        self.meshFile.write("NZONE=  %i\n" % self.n_zone)

        # Writing the full-machine mesh.
        self.writeMeshFile()
        self.meshFile.close()

    def replaceTerms(self):
        # This function replaces the boundary names of each blade row mesh to be compatible with the full-stage mesh.

        # Looping over the stages and blade rows to replace the boundary names with appropriate names.
        k = 1
        for i in range(self.n_stage):
            for j in [1, 2]:
                meshDir = self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/BladeMesh/"
                os.chdir(meshDir)
                os.system("mv mesh_out.su2 Blade_mesh.su2")

                os.system("sed -i 's|IZONE=  1"+"|IZONE= "+str(k)+"|' Blade_mesh.su2")
                os.system("sed -i 's/inflow/inflow_" + str(k) + "/g' Blade_mesh.su2")
                os.system("sed -i 's/outflow/outflow_" + str(k) + "/g' Blade_mesh.su2")
                os.system("sed -i 's/periodic1/periodic1_" + str(k) + "/g' Blade_mesh.su2")
                os.system("sed -i 's/periodic2/periodic2_" + str(k) + "/g' Blade_mesh.su2")
                k += 1

    def writeMeshFile(self):
        # This function takes the physical blade mesh from each row and appends it to the full-machine mesh so it can
        # be used for 2D blade analyses in SU2.
        for i in range(self.n_stage):
            for j in [1, 2]:
                with open(self.dir + "/Stage_"+str(i+1)+"/Bladerow_"+str(j)+"/BladeMesh/Blade_mesh.su2", "r") as Blademesh:
                    lines = Blademesh.readlines()[1:]
                    self.meshFile.writelines(lines)
                Blademesh.close()
