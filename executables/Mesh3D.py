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


class FullAnnulus:

    def __init__(self, Meangen, IN):
        self.M = Meangen
        self.model = gmsh.model
        self.factory = self.model.occ
        self.n_rows = self.M.n_stage * 2

        self.np = 3

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        self.get_annulus_geom()
        self.make_annulus()
        self.divide_sections()
        self.flowfield = self.factory.fragment([self.annulus[0]], self.shapes)
        self.annulus = self.flowfield
        self.factory.rotate(self.annulus[0], 0, 0, 0, 0, 1, 0, 0.5 * np.pi)
        self.model.occ.synchronize()
        inlet_number = self.model.getBoundary(self.annulus[0])[1]
        inlet = self.model.addPhysicalGroup(2, inlet_number)
        self.model.setPhysicalName(2, inlet, "inlet")

        outlet_number = self.model.getBoundary(self.annulus[0])[10]
        outlet = self.model.addPhysicalGroup(2, outlet_number)
        self.model.setPhysicalName(2, outlet, "outlet")

        hub_list = self.model.getBoundary(self.annulus[0])[3:13:2]
        hub = self.model.addPhysicalGroup(2, hub_list)
        self.model.setPhysicalName(2, hub, "hub")

        shroud_list = self.model.getBoundary(self.annulus[0])[0:10:2]
        shroud = self.model.addPhysicalGroup(2, shroud_list)
        self.model.setPhysicalName(2, shroud, "shroud")

        self.model.addPhysicalGroup(3, [1], 1)
        self.model.setPhysicalName(3, 1, "Flowfield")

        self.model.setColor(inlet_number, 255, 0, 0)
        self.model.setColor(outlet_number, 0, 0, 255)
        self.model.setColor(hub_list, 0, 255, 0)
        self.model.setColor(shroud_list, 127, 0, 127)

        self.model.mesh.setSize(self.model.getEntities(0), 0.03)
        self.model.mesh.setSize(self.model.getBoundary(inlet_number), 0.03)
        self.model.mesh.setSize(self.model.getBoundary(outlet_number), 0.03)

        #self.stack_sections()
        #self.factory.synchronize()
        #self.define_boundaries()
        #self.refine_mesh()

        self.model.mesh.generate(3)
        gmsh.write("fullannulus.su2")
        gmsh.fltk.run()
        gmsh.finalize()

    def refine_mesh(self):
        self.model.mesh.setSize(self.model.getEntities(0), 0.05)
        for i in range(self.n_rows):
            j = 2*i + 1
            self.model.mesh.setSize(self.model.getBoundary(self.shapes[j], False, False, True), self.refinements[j])

    def define_boundaries(self):
        annulusGroup = self.model.addPhysicalGroup(3, self.flowfield[0])
        self.model.setPhysicalName(3, annulusGroup, "Flowfield")

        inlet = self.model.addPhysicalGroup(2, [self.model.getBoundary(self.shapes[0])[1]])
        self.model.setPhysicalName(2, inlet, "inlet")

        outlet = self.model.addPhysicalGroup(2, [self.model.getBoundary(self.shapes[-1])[2]])
        self.model.setPhysicalName(2, outlet, "outlet")

        shroud_surfaces = []
        hub_surfaces = []
        shroud_surfaces.append(self.model.getBoundary(self.shapes[0])[0])
        hub_surfaces.append(self.model.getBoundary(self.shapes[0])[3])
        for i in range(1, len(self.X_hub)-1):
            print(self.shapes[i])
            section = self.model.getBoundary(self.shapes[i])
            print(section)
            shroud_surfaces.append(section[1])
            hub_surfaces.append(section[3])
        hub = self.model.addPhysicalGroup(2, hub_surfaces)
        self.model.setPhysicalName(2, hub, "hub")
        shroud = self.model.addPhysicalGroup(2, shroud_surfaces)
        self.model.setPhysicalName(2, shroud, "shroud")

    def stack_sections(self):
        self.flowfield = self.factory.fragment([self.annulus[0]], self.shapes)
        self.factory.rotate(self.flowfield[0], 0, 0, 0, 0, 1, 0, 0.5*np.pi)


    def divide_sections(self):
        self.shapes = []
        self.refinements = []
        for j in range(1, len(self.X_hub)):
            loops_hub = []
            loops_shroud = []
            for i in range(j - 1, j + 1):
                circle_hub = self.model.occ.addCircle(0, 0, self.X_hub[i], self.R_hub[i])
                loop_hub = self.model.occ.addCurveLoop([circle_hub])
                loops_hub.append(loop_hub)

                circle_shroud = self.model.occ.addCircle(0, 0, self.X_shroud[i], self.R_shroud[i])
                loop_shroud = self.model.occ.addCurveLoop([circle_shroud])
                loops_shroud.append(loop_shroud)
            self.refinements.append((self.X_hub[j] - self.X_hub[j - 1]) / self.np)
            outer_inlet = self.factory.addThruSections(loops_shroud, makeSolid=True, makeRuled=True)
            hole_inlet = self.factory.addThruSections(loops_hub, makeSolid=True, makeRuled=True)
            section = self.factory.cut(outer_inlet, hole_inlet)
            print(section)

            self.shapes.append(section[0])



    def make_annulus(self):
        loops_hub = []
        loops_shroud = []
        for i in range(len(self.X_shroud)):
            circle_hub = self.model.occ.addCircle(0, 0, self.X_hub[i], self.R_hub[i])
            loop_hub = self.model.occ.addCurveLoop([circle_hub])
            loops_hub.append(loop_hub)

            circle_shroud = self.model.occ.addCircle(0, 0, self.X_shroud[i], self.R_shroud[i])
            loop_shroud = self.model.occ.addCurveLoop([circle_shroud])
            loops_shroud.append(loop_shroud)

        outer_annulus = self.factory.addThruSections(loops_shroud, makeSolid=True, makeRuled=True)
        inner_annulus = self.factory.addThruSections(loops_hub, makeSolid=True, makeRuled=True)

        self.annulus = self.factory.cut(outer_annulus, inner_annulus)


    def get_annulus_geom(self):
        X_LE = self.M.X_LE
        R_LE = self.M.Z_LE
        X_TE = self.M.X_TE
        R_TE = self.M.Z_TE

        # Calculating number of rows.
        n_rows = len(X_LE[0, :])

        # Creating a list of chord lengths which will be used for local cell size calculation.
        chords = np.zeros(2 * self.M.n_stage)
        for i in range(self.M.n_stage):
            if self.M.machineType == 'C':
                chords[2 * i] += self.M.chord_R[i]
                chords[2 * i + 1] += self.M.chord_S[i]
            else:
                chords[2 * i] += self.M.chord_S[i]
                chords[2 * i + 1] += self.M.chord_R[i]
        self.Chords = chords
        X_inlet = min(X_LE[:, 0]) - 2 * chords[0]
        X_outlet = max(X_TE[:, -1]) + 2*chords[-1]

        self.X_shroud = [X_inlet]
        self.R_shroud = [R_LE[-1, 0]]
        self.X_hub = [X_inlet]
        self.R_hub = [R_LE[0, 0]]
        for i in range(n_rows):
            self.X_shroud.append(X_LE[-1, i])
            self.X_shroud.append(X_TE[-1, i])
            self.X_hub.append(X_LE[0, i])
            self.X_hub.append(X_TE[0, i])

            self.R_shroud.append(R_LE[-1, i])
            self.R_hub.append(R_LE[0, i])
            self.R_shroud.append(R_TE[-1, i])
            self.R_hub.append(R_TE[0, i])
        self.X_shroud.append(X_outlet)
        self.X_hub.append(X_outlet)
        self.R_shroud.append(R_TE[-1, -1])
        self.R_hub.append(R_TE[0, -1])

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
        chords = np.zeros(2 * self.M.n_stage)
        for i in range(self.M.n_stage):
            if self.M.machineType == 'C':
                chords[2 * i] += self.M.chord_R[i]
                chords[2 * i + 1] += self.M.chord_S[i]
            else:
                chords[2 * i] += self.M.chord_S[i]
                chords[2 * i + 1] += self.M.chord_R[i]

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
        X_inlet = min(X_LE[:, 0] - 2 * chords[0])
        X_p.append(X_inlet)
        # Calculating y-coordinate.
        Y_p.append(R_LE[0, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        # Calculating z-coordinate.
        Z_p.append(R_LE[0, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        # Inlet cell size is defined as the first blade row cell size multiplied by the inlet cell size factor.
        # L_c.append(self.inlet_fac * chords[0] / self.n_point)
        L_c.append(self.inlet_fac * chords[0] / self.n_point)
        i_point += 1

        # The second point is located at the inlet at the shroud(above the first point).
        X_p.append(X_inlet)
        Y_p.append(R_LE[-1, 0] * np.sin(0.5 * self.wedge * np.pi / 180))
        Z_p.append(R_LE[-1, 0] * np.cos(0.5 * self.wedge * np.pi / 180))
        points.append(i_point)
        # L_c.append(self.inlet_fac * chords[0] / self.n_point)
        L_c.append(self.inlet_fac * chords[0] / self.n_point)
        i_point += 1

        # Looping over the blade rows to create the points on the hub and shroud surface.
        for j in range(n_rows):
            # Defining point at the leading edge of the blade row at the hub.
            X_p.append(X_LE[0, j])
            Y_p.append(R_LE[0, j] * np.sin(0.5 * self.wedge * np.pi/180))
            Z_p.append(R_LE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append(chords[j] / self.n_point)
            i_point += 1

            # Defining point at the leading edge of the blade row at the shroud.
            X_p.append(X_LE[-1, j])
            Y_p.append(R_LE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
            Z_p.append(R_LE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append(chords[j] / self.n_point)
            i_point += 1

            # Defining point at the trailing edge of the blade row at the hub.
            X_p.append(X_TE[0, j])
            Y_p.append(R_TE[0, j] * np.sin(0.5 * self.wedge * np.pi / 180))
            Z_p.append(R_TE[0, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append(chords[j] / self.n_point)
            i_point += 1

            # Defining point at the trailing edge of the blade row at the shroud.
            X_p.append(X_TE[-1, j])
            Y_p.append(R_TE[-1, j] * np.sin(0.5 * self.wedge * np.pi / 180))
            Z_p.append(R_TE[-1, j] * np.cos(0.5 * self.wedge * np.pi / 180))
            points.append(i_point)
            L_c.append(chords[j] / self.n_point)
            i_point += 1

        # The outlet is located 2 axial chords of the last blade row downstream of its trailing edge.
        X_outlet = max(X_TE[:, -1] + 2 * chords[-1])
        X_p.append(X_outlet)
        Y_p.append(Y_p[-2])
        Z_p.append(Z_p[-2])
        points.append(i_point)
        # L_c.append(self.outlet_fac * chords[-1] / self.n_point)
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