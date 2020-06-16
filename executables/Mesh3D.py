# --------------------------------------------------------------------------------------------------------------- #
# This program contains the classes that construct the 3D mesh for BFM analysis and in the future, physical blade
# analysis. Gmesh is used for the 3D BFM mesh generation and should be added to the user's PATH and PYTHOPATH
# in the right way. For periodic mesh generation, SU2_PERIO should be installed and added to the PATH according to
# installation instructions. The gmsh module should be installed for the user's version of python. 3D physical blade
# mesh generation is not yet implemented.
# --------------------------------------------------------------------------------------------------------------- #

import gmsh
import numpy as np
import os

# Class capable of generating a 3D mesh suitable for BFM analysis in SU2.
class Gmesh3D:
    wedge =         4    # 3D wedge angle in degrees. Should be lower than 180.
    n_sec =         4    # Number of sections in tangential direction in the wedge.
    n_point =       20   # Number of points in axial direction for each blade row.
    inlet_fac =     1.0  # Inlet cell size factor. High numbers mean a more coarse inlet grid.
    outlet_fac =    1.0  # Outlet cell size factor.
    fileName =      "3DBFM.su2"  # Mesh file name.
    pointID =       None  # List of point identification numbers.
    xCoords =       None  # List of x-coordinates of the mesh points.
    yCoords =       None  # List of y-coordinates of the mesh points.
    zCoords =       None  # List of z-coordinates of the mesh points.
    lC =            None  # List of mesh cell sizes.

    lineID =        None  # List of line identification numbers.
    lineStart =     None  # List of line start point identification numbers.
    lineEnd =       None  # List of line end point identification numbers.

    Rev =           None  # Revolution class.

    def __init__(self, Meangen, IN):
        # Storing Meangen class.
        self.M = Meangen

        # Reading mesh parameters from input file.
        self.wedge = IN["WEDGE"][0]
        self.n_sec = int(IN["SECTIONS_PER_DEGREE"][0]) * self.wedge
        self.n_point = int(IN["AXIAL_POINTS"][0])
        self.inlet_fac = IN["INLET_FACTOR"][0]
        self.outlet_fac = IN["OUTLET_FACTOR"][0]

        # Initiating Gmesh.
        self.model = gmsh.model
        self.factory = self.model.geo
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        self.model.add("3DBFM")

        # Creating all the mesh points.
        self.makePoints()

        # Connecting the mesh points to create lines.
        self.makeLines()

        # Making periodic plane.
        self.makePlane()

        # Revolving plane around rotation axis to create 3D wedge.
        self.revolve()

        # Setting names to boundaries.
        self.nameBoundaries()

        # Creating 3D mesh.
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.write(self.fileName)

        # In case of selection of mesh visualization option, the gmesh GUI will display the 3D BFM mesh.
        if IN["PLOT_MESH"] == 'YES':
            gmsh.fltk.run()

        # Transforming mesh to periodic mesh.
        print("Building periodic mesh...")
        self.makePerio()
        print("Done!")

    def makePerio(self):
        # This function modifies the mesh to be periodic through SU2_PERIO.
        # Creating SU2_PERIO configuration file.
        file = open("createPerio.cfg", "w+")

        # Writing periodic boundary command and specifying wedge angle.
        file.write("MARKER_PERIODIC= (periodic_1, periodic_2, 0.0, 0.0, 0.0, " + str(
            float(self.wedge)) + ", 0.0, 0.0, 0.0, 0.0, 0.0)\n")
        file.write("MESH_FILENAME= "+self.fileName+"\n")
        file.write("MESH_FORMAT= SU2\n")
        file.write("MESH_OUT_FILENAME= "+self.fileName[:-4]+"_perio.su2\n")
        file.close()

        # Executing SU2_PERIO to create periodic mesh and storing output in output file.
        os.system("SU2_PERIO < createPerio.cfg > SU2_PERIO.out")

    def nameBoundaries(self):
        # This function gives names to all the boundaries of the 3D mesh so boundary conditions can be assigned in the
        # SU2 configuration file.

        # Specifying periodic boundaries.
        # First periodic boundary is the plane created by makePlane().
        periodic_1 = self.model.addPhysicalGroup(2, [1])
        self.model.setPhysicalName(2, periodic_1, "periodic_1")
        # Second periodic boundary is the first plane in the revolution volume list.
        periodic_2 = self.model.addPhysicalGroup(2, [self.Rev[0][1]])
        self.model.setPhysicalName(2, periodic_2, "periodic_2")

        # Setting the revolution volume as a 3D physical group.
        self.model.addPhysicalGroup(3, [1], 1)
        self.model.setPhysicalName(3, 1, "FlowField")

        # Going over the other surfaces of the revolution volume and naming them accordingly.
        i = 2

        # Naming the inlet boundary.
        inlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
        self.model.setPhysicalName(2, inlet, "inlet")
        i += 1

        # Appending all planes on the shroud side of the volume to a single shroud boundary.
        shroud_list = []
        for j in range(i, i + 3 * self.M.n_stage + 2):
            shroud_list.append(self.Rev[j][1])
            i += 1
        shroud = self.model.addPhysicalGroup(2, shroud_list)
        self.model.setPhysicalName(2, shroud, "shroud")

        # Naming the outlet boundary
        outlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
        self.model.setPhysicalName(2, outlet, "outlet")
        i += 1

        # Appending all planes on the shroud side of the volume to a single shroud boundary.
        hub_list = []
        for j in range(i, i + 3 * self.M.n_stage + 2):
            hub_list.append(self.Rev[j][1])
            i += 1
        hub = self.model.addPhysicalGroup(2, hub_list)
        self.model.setPhysicalName(2, hub, "hub")

    def revolve(self):
        # This function revolves the 2D periodic plane around the rotation axis to create a wedge.
        # The center point of the revolution is located at the origin and the rotation axis is set to be the x-axis.
        self.Rev = self.factory.revolve([(2, 1)], 0, 0, 0, 1, 0, 0, self.wedge*np.pi/180, [self.n_sec])

    def makePlane(self):
        # This function creates the first periodic plane by adding all the points, creating the lines and appending
        # them to a curveloop to create a 2D plane mesh.

        # Adding all the points to the mesh.
        for i in range(len(self.pointID)):
            self.factory.addPoint(self.xCoords[i], self.yCoords[i], self.zCoords[i], self.lC[i], self.pointID[i])

        # Connecting the points by creating lines.
        for i in range(len(self.lineID)):
            self.factory.addLine(self.lineStart[i], self.lineEnd[i], self.lineID[i])

        # Initiating curveloop and specifying surface.
        self.factory.addCurveLoop(self.lineID, 1)
        self.factory.addPlaneSurface([1], 1)

    def makeLines(self):
        # This function sets up the start and end points of the lines that make up the periodic plane.

        # Defining property lists.
        i_line = 1       # Current line ID
        start_line = []  # List containing line start points.
        end_line = []    # List containing line end points.
        lines = []       # List containing line tags.

        # Looping over points in clockwise direction, starting at the hub at the inlet.
        # The inlet is set as the first line.
        start_line.append(self.pointID[0])
        end_line.append(self.pointID[1])
        lines.append(i_line)
        i_line += 1

        # Looping over the points on the shroud line.
        for i in range(1, len(self.pointID)-1, 2):
            start_line.append(self.pointID[i])
            end_line.append(self.pointID[i+2])
            lines.append(i_line)
            i_line += 1

        # Going down the outlet line.
        start_line.append(self.pointID[-1])
        end_line.append(self.pointID[-2])
        lines.append(i_line)
        i_line += 1

        # Moving back to the starting point via the hub line.
        for i in range(len(self.pointID)-2, 0, -2):
            start_line.append(self.pointID[i])
            end_line.append(self.pointID[i-2])
            lines.append(i_line)
            i_line += 1

        # Storing line tags, starting points and end points.
        self.lineID = lines
        self.lineStart = start_line
        self.lineEnd = end_line

    def makePoints(self):
        # This function sets up the points which make up the periodic plane of the BMF mesh.

        # Extracting the leading and trailing edge coordinates from Meangen.
        X_LE = self.M.X_LE
        R_LE = self.M.Z_LE
        X_TE = self.M.X_TE
        R_TE = self.M.Z_TE

        # Calculating number of rows.
        n_rows = len(X_LE[0, :])

        # Creating a list of chord lengths which will be used for local cell size calculation.
        chords = []
        for i in range(self.M.n_stage):
            if self.M.machineType == 'C':
                chords.append(self.M.chord_R[i])
                chords.append(self.M.chord_S[i])
            else:
                chords.append(self.M.chord_S[i])
                chords.append(self.M.chord_R[i])

        # Defining lists for x, y and z coordinates of the points, the cell size for each point and point tags.
        X_p = []
        Y_p = []
        Z_p = []
        L_c = []
        points = []

        # The first point is located at the domain inlet at the hub, which is 2 axial chords in front of the first
        # blade row leading edge.
        i_point = 1
        # Calculating x-coordinate.
        X_inlet = min(X_LE[0, 0] - 2 * chords[0], X_LE[-1, 0] - 2 * chords[0])
        X_p.append(X_inlet)
        # Calculating y-coordinate.
        Y_p.append(R_LE[0, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        # Calculating z-coordinate.
        Z_p.append(R_LE[0, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        # Inlet cell size is defined as the first blade row cell size multiplied by the inlet cell size factor.
        L_c.append(self.inlet_fac * chords[0] / self.n_point)
        i_point += 1

        # The second point is located at the inlet at the shroud(above the first point).
        X_p.append(X_inlet)
        Y_p.append(R_LE[-1, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        Z_p.append(R_LE[-1, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        L_c.append(self.inlet_fac * chords[0] / self.n_point)
        i_point += 1

        # Looping over the blade rows to create the points on the hub and shroud surface.
        for j in range(n_rows):
            # Defining point at the leading edge of the blade row at the hub.
            X_p.append(X_LE[0, j])
            Y_p.append(R_LE[0, j] * np.sin(0.5 * self.wedge * np.pi/180))
            Z_p.append(R_LE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append((X_TE[0, j] - X_LE[0, j]) / self.n_point)
            i_point += 1

            # Defining point at the leading edge of the blade row at the shroud.
            X_p.append(X_LE[-1, j])
            Y_p.append(R_LE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
            Z_p.append(R_LE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append((X_TE[-1, j] - X_LE[-1, j]) / self.n_point)
            i_point += 1

            # Defining point at the trailing edge of the blade row at the hub.
            X_p.append(X_TE[0, j])
            Y_p.append(R_TE[0, j] * np.sin(0.5 * self.wedge * np.pi / 180))
            Z_p.append(R_TE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append((X_TE[0, j] - X_LE[0, j]) / self.n_point)
            i_point += 1

            # Defining point at the trailing edge of the blade row at the shroud.
            X_p.append(X_TE[-1, j])
            Y_p.append(R_TE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
            Z_p.append(R_TE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append((X_TE[-1, j] - X_LE[-1, j]) / self.n_point)
            i_point += 1

        # The outlet is located 2 axial chords of the last blade row downstream of its trailing edge.
        X_outlet = max(X_p[-2] + 2 * chords[-1], X_p[-2] + 2 * chords[-1])
        X_p.append(X_outlet)
        Y_p.append(Y_p[-2])
        Z_p.append(Z_p[-2])
        points.append(i_point)
        L_c.append(self.outlet_fac * chords[-1] / self.n_point)
        i_point += 1
        X_p.append(X_outlet)
        Y_p.append(Y_p[-2])
        Z_p.append(Z_p[-2])
        points.append(i_point)
        L_c.append(self.outlet_fac * chords[-1] / self.n_point)
        i_point += 1

        # Storing point tags and coordinates.
        self.pointID = points
        self.xCoords = X_p
        self.yCoords = Y_p
        self.zCoords = Z_p
        self.lC = L_c