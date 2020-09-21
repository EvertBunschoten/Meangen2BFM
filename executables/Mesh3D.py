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
        #self.divide_sections()
        #self.flowfield = self.factory.fragment([self.annulus[0]], self.shapes)
        #self.annulus = self.flowfield
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
        #self.factory.rotate(self.flowfield[0], 0, 0, 0, 0, 1, 0, 0.5*np.pi)


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

# Class capable of generating a 3D mesh suitable for axisymmetric BFM analysis in SU2.
class Gmesh3D:
    wedge =         1    # 3D wedge angle in degrees. Should be lower than 180.
    n_sec =         1    # Number of sections in tangential direction in the wedge.
    n_point =       20   # Number of points in axial direction for each blade row.
    inlet_fac =     1.0  # Inlet cell size factor. High numbers mean a more coarse inlet grid.
    outlet_fac =    1.0  # Outlet cell size factor.
    fileName =      "3DBFM.su2"  # Mesh file name.
    pointID =       None  # List of point identification numbers.
    xCoords =       None  # List of x-coordinates of the mesh points.
    yCoords =       None  # List of y-coordinates of the mesh points.
    zCoords =       None  # List of z-coordinates of the mesh points.
    lC =            None  # List of mesh cell sizes.
    rot_axis =      [0, 0, 1]
    lineID =        None  # List of line identification numbers.
    lineStart =     None  # List of line start point identification numbers.
    lineEnd =       None  # List of line end point identification numbers.

    Rev =           None  # Revolution class.


    def __init__(self, Meangen, IN):
        # Storing Meangen class.
        self.M = Meangen

        # Reading mesh parameters from input file.

        self.wedge = IN["WEDGE"][0]     # Importing mesh wedge angle
        self.n_sec = int(IN["SECTIONS_PER_DEGREE"][0]) * self.wedge     # Calculating number of tangential nodes
        self.n_point = int(IN["AXIAL_POINTS"][0])       # Importing axial node count
        self.inlet_fac = IN["INLET_FACTOR"][0]          # Importing inlet cell size factor
        self.outlet_fac = IN["OUTLET_FACTOR"][0]        # Importing outlet cell size factor
        self.rot_axis = IN["Rotation_axis"]

        self.Coords = type('', (), {})()
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

        # Applying mesh refinement along the lines in the mesh
        self.refineLines()

        # Generate the 3D mesh geometry
        gmsh.model.mesh.generate(3)

        # Saving mesh as su2 file
        gmsh.write(self.fileName)

        # Transforming mesh to periodic mesh.
        print("Building periodic mesh...")
        self.makePerio()
        print("Done!")
        # In case of selection of mesh visualization option, the gmesh GUI will display the 3D BFM mesh.
        if IN["PLOT_MESH"] == 'YES':
            self.plotmesh()
        gmsh.finalize()
    def plotmesh(self):
        v = gmsh.view.add("comments")
        inlet_coords = np.sum(self.Coords.Coords_inlet, axis=0)/np.shape(self.Coords.Coords_inlet)[0]
        outlet_coords = np.sum(self.Coords.Coords_outlet, axis=0) / np.shape(self.Coords.Coords_outlet)[0]
        hub_coords = np.sum(self.Coords.Coords_hub, axis=0) / np.shape(self.Coords.Coords_hub)[0]
        shroud_coords = np.sum(self.Coords.Coords_shroud, axis=0) / np.shape(self.Coords.Coords_shroud)[0]
        perio_1_coords = np.sum(self.Coords.Coords_periodic_1, axis=0) / np.shape(self.Coords.Coords_periodic_1)[0]
        perio_2_coords = np.sum(self.Coords.Coords_periodic_2, axis=0) / np.shape(self.Coords.Coords_periodic_2)[0]

        gmsh.view.addListDataString(v, [i for i in inlet_coords], ["inlet"], ["Align", "Center", "Font", "Helvetica"])
        gmsh.view.addListDataString(v, [i for i in outlet_coords], ["outlet"], ["Align", "Center", "Font", "Helvetica"])
        gmsh.view.addListDataString(v, [i for i in hub_coords], ["hub"], ["Align", "Center", "Font", "Helvetica"])
        gmsh.view.addListDataString(v, [i for i in shroud_coords], ["shroud"], ["Align", "Center", "Font", "Helvetica"])
        gmsh.view.addListDataString(v, [i for i in perio_2_coords], ["periodic_2"],
                                    ["Align", "Center", "Font", "Helvetica"])
        gmsh.view.addListDataString(v, [i for i in perio_1_coords], ["periodic_1"], ["Align", "Center", "Font", "Helvetica"])

        gmsh.fltk.run()
    def makePerio(self):
        # This function modifies the mesh to be periodic through SU2_PERIO.
        # Creating SU2_PERIO configuration file.
        file = open("createPerio.cfg", "w+")

        # Writing periodic boundary command and specifying wedge angle.
        file.write("MARKER_PERIODIC= (periodic_1, periodic_2, 0.0, 0.0, 0.0, 0.0, 0.0, " + str(
            float(self.wedge)) + ", 0.0, 0.0, 0.0)\n")
        file.write("MESH_FILENAME= "+self.fileName+"\n")
        file.write("MESH_FORMAT= SU2\n")
        file.write("MESH_OUT_FILENAME= "+self.fileName[:-4]+"_perio.su2\n")
        file.close()

        # Executing SU2_PERIO to create periodic mesh and storing output in output file.
        os.system("SU2_PERIO createPerio.cfg")

    def nameBoundaries(self):
        # This function gives names to all the boundaries of the 3D mesh so boundary conditions can be assigned in the
        # SU2 configuration file.

        # Specifying periodic boundaries.
        # First periodic boundary is the plane created by makePlane().
        periodic_1 = self.model.addPhysicalGroup(2, [1])
        self.model.setPhysicalName(2, periodic_1, "periodic_1")
        self.model.setColor([(2, 1)], 255, 255, 0)
        # Second periodic boundary is the first plane in the revolution volume list.
        periodic_2 = self.model.addPhysicalGroup(2, [self.Rev[0][1]])
        self.model.setPhysicalName(2, periodic_2, "periodic_2")
        self.model.setColor([self.Rev[-1]], 255, 102, 0)

        # Setting the revolution volume as a 3D physical group.
        self.model.addPhysicalGroup(3, [1], 1)
        self.model.setPhysicalName(3, 1, "FlowField")
        self.model.setColor((3, 1), 0, 0, 0)
        # Going over the other surfaces of the revolution volume and naming them accordingly.
        i = 2

        # Naming the inlet boundary.
        inlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
        self.model.setPhysicalName(2, inlet, "inlet")
        self.model.setColor([self.Rev[i]], 255, 0, 0)
        i += 1

        # Appending all planes on the shroud side of the volume to a single shroud boundary.
        shroud_list = []
        for j in range(i, i + 3 * self.M.n_stage + 2):
            shroud_list.append(self.Rev[j][1])
            i += 1
        shroud = self.model.addPhysicalGroup(2, shroud_list)
        self.model.setPhysicalName(2, shroud, "shroud")
        self.model.setColor([(2, j) for j in shroud_list], 0, 255, 0)
        # Naming the outlet boundary
        outlet = self.model.addPhysicalGroup(2, [self.Rev[i][1]])
        self.model.setPhysicalName(2, outlet, "outlet")
        self.model.setColor([self.Rev[i]], 0, 0, 255)
        i += 1

        # Appending all planes on the shroud side of the volume to a single shroud boundary.
        hub_list = []
        for j in range(i, i + 3 * self.M.n_stage + 2):
            hub_list.append(self.Rev[j][1])
            i += 1
        hub = self.model.addPhysicalGroup(2, hub_list)
        self.model.setPhysicalName(2, hub, "hub")
        self.model.setColor([(2, j) for j in hub_list], 255, 0, 255)
    def refineLines(self):
        # This function sets the number of nodes along each of the curves bounding the geometry. This allows for
        # refinements within the bladed regions along the leading and trailing edges of the blades, while the mesh is
        # kept relatively coarse near the inlet and outlet patches

        # The application of the mesh refinement process follows the natural progression of the machine, starting with
        # the inlet farfield.
        # Setting node count along hub and shroud boundaries.
        self.model.mesh.setTransfiniteCurve(self.lines_hub[0], self.n_point)
        self.model.mesh.setTransfiniteCurve(self.lines_shroud[0], self.n_point)
        #
        for i in range(1, len(self.lines_hub)-1):
            if i % 2 == 0:
                self.model.mesh.setTransfiniteCurve(self.lines_hub[i], max([int(self.n_point*self.length_hub[i]/self.length_hub[i-1]), 3]))
                self.model.mesh.setTransfiniteCurve(self.lines_shroud[i], max([int(self.n_point*self.length_shroud[i]/self.length_shroud[i-1]), 3]))
            else:
                self.model.mesh.setTransfiniteCurve(self.lines_hub[i], self.n_point)
                self.model.mesh.setTransfiniteCurve(self.lines_shroud[i], self.n_point)
        self.model.mesh.setTransfiniteCurve(self.lines_hub[-1], self.n_point)
        self.model.mesh.setTransfiniteCurve(self.lines_shroud[-1], self.n_point)

        for i in range(len(self.lines_rad)):
            self.model.mesh.embed(1, [self.lines_rad[i]], 2, 1)
            if i == 0:
                self.model.mesh.setTransfiniteCurve(self.lines_rad[i], int(self.n_point/self.inlet_fac))
            elif i == len(self.lines_rad)-1:
                self.model.mesh.setTransfiniteCurve(self.lines_rad[i], int(self.n_point / self.outlet_fac))
            else:
                if i % 2 == 0:
                    L_av = 0.5 * (self.length_hub[i-1] + self.length_shroud[i-1])
                else:
                    L_av = 0.5 * (self.length_hub[i] + self.length_shroud[i])
                self.model.mesh.setTransfiniteCurve(self.lines_rad[i], max([int(self.n_point*self.length_rad[i]/L_av), 3]))
    def revolve(self):
        # This function revolves the 2D periodic plane around the rotation axis to create a wedge.
        # The center point of the revolution is located at the origin and the rotation axis is set to be the x-axis.
        self.Rev = self.factory.revolve([(2, 1)], 0, 0, 0, 0, 0, 1, self.wedge*np.pi/180, [self.n_sec])

    def makePlane(self):
        loop = []
        loop.append(-self.lines_rad[0])
        for i in self.lines_hub:
            loop.append(i)
        loop.append(self.lines_rad[-1])
        for i in self.lines_shroud[::-1]:
            loop.append(-i)
        self.factory.addCurveLoop(loop, 1)
        self.factory.addPlaneSurface([1], 1)
    def makeLines(self):
        lines_hub = []
        length_hub = []
        lines_shroud = []
        length_shroud = []
        lines_rad = []
        length_rad = []
        #print(np.sqrt(np.sum(np.array(self.coords_hub[:, 1] - self.coords_hub[:, 0])**2)))
        i_line = 1
        for i in range(len(self.points_hub)-1):
            self.factory.addLine(self.points_hub[i], self.points_hub[i+1], i_line)
            lines_hub.append(i_line)
            length_hub.append(np.sqrt(np.sum(np.array(self.coords_hub[:, i+1] - self.coords_hub[:, i])**2)))
            i_line += 1
        for i in range(len(self.points_shroud)-1):
            self.factory.addLine(self.points_shroud[i], self.points_shroud[i+1], i_line)
            lines_shroud.append(i_line)
            length_shroud.append(np.sqrt(np.sum(np.array(self.coords_shroud[:, i + 1] - self.coords_shroud[:, i]) ** 2)))
            i_line += 1
        for i in range(len(self.points_hub)):
            self.factory.addLine(self.points_hub[i], self.points_shroud[i], i_line)
            lines_rad.append(i_line)
            length_rad.append(np.sqrt(np.sum(np.array(self.coords_shroud[:, i] - self.coords_hub[:, i]) ** 2)))
            i_line += 1
        self.lines_hub = lines_hub
        self.lines_shroud = lines_shroud
        self.lines_rad = lines_rad

        self.length_hub = length_hub
        self.length_shroud = length_shroud
        self.length_rad = length_rad
    def makePoints(self):
        # Extracting the leading and trailing edge coordinates from Meangen.
        X_LE = self.M.X_LE
        R_LE = self.M.Z_LE
        X_TE = self.M.X_TE
        R_TE = self.M.Z_TE

        # Calculating number of rows.
        n_rows = len(X_LE[0, :])

        x_min = min(X_LE[:, 0] - 2*(X_TE[:, 0] - X_LE[:, 0]))
        x_max = max(X_TE[:, -1] + 2*(X_TE[:, -1] - X_LE[:, -1]))

        i_point = 1
        Z_hub = []
        X_hub = []
        Y_hub = []
        Z_shroud = []
        X_shroud = []
        Y_shroud = []

        points_hub = []
        points_shroud = []

        theta = 0.5*self.wedge*np.pi/180
        Z_hub.append(x_min)
        X_hub.append(R_LE[0, 0]*np.sin(theta))
        Y_hub.append(R_LE[0, 0] * np.cos(theta))
        points_hub.append(i_point)
        i_point += 1

        for i in range(n_rows):
            Z_hub.append(X_LE[0, i])
            X_hub.append(R_LE[0, i] * np.sin(theta))
            Y_hub.append(R_LE[0, i] * np.cos(theta))
            points_hub.append(i_point)
            i_point += 1

            Z_hub.append(X_TE[0, i])
            X_hub.append(R_TE[0, i] * np.sin(theta))
            Y_hub.append(R_TE[0, i] * np.cos(theta))
            points_hub.append(i_point)
            i_point += 1

        Z_hub.append(x_max)
        X_hub.append(R_TE[0, -1] * np.sin(theta))
        Y_hub.append(R_TE[0, -1] * np.cos(theta))
        points_hub.append(i_point)
        i_point += 1

        Z_shroud.append(x_min)
        X_shroud.append(R_LE[-1, 0] * np.sin(theta))
        Y_shroud.append(R_LE[-1, 0] * np.cos(theta))
        points_shroud.append(i_point)
        i_point += 1
        for i in range(n_rows):
            Z_shroud.append(X_LE[-1, i])
            X_shroud.append(R_LE[-1, i] * np.sin(theta))
            Y_shroud.append(R_LE[-1, i] * np.cos(theta))
            points_shroud.append(i_point)
            i_point += 1

            Z_shroud.append(X_TE[-1, i])
            X_shroud.append(R_TE[-1, i] * np.sin(theta))
            Y_shroud.append(R_TE[-1, i] * np.cos(theta))
            points_shroud.append(i_point)
            i_point += 1
        Z_shroud.append(x_max)
        X_shroud.append(R_TE[-1, -1] * np.sin(theta))
        Y_shroud.append(R_TE[-1, -1] * np.cos(theta))
        points_shroud.append(i_point)
        i_point += 1

        self.Coords.Coords_inlet = np.array([[X_hub[0], Y_hub[0], Z_hub[0]],
                                           [X_shroud[0], Y_shroud[0], Z_shroud[0]],
                                           [-X_hub[0], Y_hub[0], Z_hub[0]],
                                           [-X_shroud[0], Y_shroud[0], Z_shroud[0]]])
        self.Coords.Coords_outlet = np.array([[X_hub[-1], Y_hub[-1], Z_hub[-1]],
                                           [X_shroud[-1], Y_shroud[-1], Z_shroud[-1]],
                                           [-X_hub[-1], Y_hub[-1], Z_hub[-1]],
                                           [-X_shroud[-1], Y_shroud[-1], Z_shroud[-1]]])
        self.Coords.Coords_hub = np.transpose(np.array([X_hub + [-x for x in X_hub], Y_hub + Y_hub, Z_hub + Z_hub]))
        self.Coords.Coords_shroud = np.transpose(np.array([X_shroud + [-x for x in X_shroud], Y_shroud + Y_shroud, Z_shroud + Z_shroud]))
        self.Coords.Coords_periodic_1 = np.transpose(np.array([X_hub + X_shroud, Y_hub + Y_shroud, Z_hub + Z_shroud]))
        self.Coords.Coords_periodic_2 = np.transpose(np.array([[-x for x in X_hub] + [-x for x in X_shroud],
                                                               [x for x in Y_hub] + [x for x in Y_shroud],
                                                               [x for x in Z_hub] + [x for x in Z_shroud]]))

        self.points_hub = points_hub
        self.coords_hub = np.mat([X_hub, Y_hub, Z_hub])
        self.points_shroud = points_shroud
        self.coords_shroud = np.mat([X_shroud, Y_shroud, Z_shroud])
        for i in range(len(X_hub)):
            self.factory.addPoint(X_hub[i], Y_hub[i], Z_hub[i], 0.01, points_hub[i])
            self.factory.addPoint(X_shroud[i], Y_shroud[i], Z_shroud[i], 0.01, points_shroud[i])

class Gmesh2D:
    wedge = 1.0
    M = None
    Np = 10
    points_perio_1 = None
    points_perio_2 = None
    lines_perio_1 = None
    lines_perio_2 = None
    lines_rad = None
    fileName = "2DBFM"

    def __init__(self, Meangen, IN):
        self.wedge = IN["WEDGE"][0]
        self.Np = int(IN["AXIAL_POINTS"][0])
        self.M = Meangen

        self.model = gmsh.model
        self.factory = self.model.geo
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        self.model.add("2DBFM")

        self.makepoints()
        self.makelines()
        self.makesurface()
        self.factory.synchronize()
        self.refinelines()
        gmsh.model.mesh.generate(3)
        self.nameboundaries()



        gmsh.write(self.fileName+".su2")
        self.writeperio()
        if IN["PLOT_MESH"] == 'YES':
            gmsh.fltk.run()

    def writeperio(self):

        file = open("createPerio.cfg", "w+")

        # Writing periodic boundary command and specifying wedge angle.
        file.write("MARKER_PERIODIC= (periodic_1, periodic_2, 0.0, 0.0, 0.0, 0.0, 0.0, " + str(
            self.wedge) + ", 0.0, 0.0, 0.0)\n")
        file.write("MESH_FILENAME= " + self.fileName + ".su2\n")
        file.write("MESH_FORMAT= SU2\n")
        file.write("MESH_OUT_FILENAME= " + self.fileName + "_perio.su2\n")
        file.close()

        # Executing SU2_PERIO to create periodic mesh and storing output in output file.
        os.system("SU2_PERIO createPerio.cfg > SU2Output/SU2_perio.out")
    def nameboundaries(self):
        self.model.addPhysicalGroup(2, [1], 1)
        self.model.setPhysicalName(2, 1, "Flowfield")
        periodic_1 = self.model.addPhysicalGroup(1, self.lines_perio_1)
        self.model.setPhysicalName(1, periodic_1, "periodic_1")
        periodic_2 = self.model.addPhysicalGroup(1, self.lines_perio_2)
        self.model.setPhysicalName(1, periodic_2, "periodic_2")

        inlet = self.model.addPhysicalGroup(1, [self.lines_rad[0]])
        self.model.setPhysicalName(1, inlet, "inlet")

        outlet = self.model.addPhysicalGroup(1, [self.lines_rad[-1]])
        self.model.setPhysicalName(1, outlet, "outlet")

    def refinelines(self):
        self.model.mesh.setTransfiniteCurve(self.lines_perio_1[0], self.Np)
        self.model.mesh.setTransfiniteCurve(self.lines_perio_2[0], self.Np)
        for i in range(1, len(self.lines_perio_1) - 1):
            if i % 2 == 0:
                self.model.mesh.setTransfiniteCurve(self.lines_perio_1[i], int(self.Np / 2))
                self.model.mesh.setTransfiniteCurve(self.lines_perio_2[i], int(self.Np / 2))
            else:
                self.model.mesh.setTransfiniteCurve(self.lines_perio_1[i], self.Np)
                self.model.mesh.setTransfiniteCurve(self.lines_perio_2[i], self.Np)
        self.model.mesh.setTransfiniteCurve(self.lines_perio_1[-1], self.Np)
        self.model.mesh.setTransfiniteCurve(self.lines_perio_2[-1], self.Np)

    def makesurface(self):
        loop = []
        loop.append(self.lines_rad[0])
        loop += [i for i in self.lines_perio_2]
        loop.append(-self.lines_rad[-1])
        loop += [-i for i in self.lines_perio_1]

        self.factory.addCurveLoop(loop, 1)
        self.factory.addPlaneSurface([1], 1)

    def makelines(self):
        lines_perio_1 = []
        lines_perio_2 = []
        lines_rad = []
        i_line = 1

        for i in range(len(self.points_perio_1) - 1):
            self.factory.addLine(self.points_perio_1[i], self.points_perio_1[i + 1], i_line)
            lines_perio_1.append(i_line)
            i_line += 1

            self.factory.addLine(self.points_perio_2[i], self.points_perio_2[i + 1], i_line)
            lines_perio_2.append(i_line)
            i_line += 1

            self.factory.addLine(self.points_perio_1[i], self.points_perio_2[i], i_line)
            lines_rad.append(i_line)
            i_line += 1

        self.factory.addLine(self.points_perio_1[-1], self.points_perio_2[-1], i_line)
        lines_rad.append(i_line)

        self.lines_perio_1 = lines_perio_1
        self.lines_perio_2 = lines_perio_2
        self.lines_rad = lines_rad

    def makepoints(self):
        X_LE = self.M.X_LE
        X_TE = self.M.X_TE
        R_LE = self.M.Z_LE
        R_TE = self.M.Z_TE

        x_inlet = X_LE[1, 0] - 2 * (X_TE[1, 0] - X_LE[1, 0])
        x_outlet = X_TE[1, -1] + 2 * (X_TE[1, -1] - X_LE[1, -1])

        i_point = 1
        points_perio_1 = []
        points_perio_2 = []

        theta = self.wedge*np.pi/180

        self.factory.addPoint(R_LE[1, 0] * np.sin(0.5 * theta), R_LE[1, 0] * np.cos(0.5 * theta), x_inlet, 1.0, i_point)
        points_perio_1.append(i_point)
        i_point += 1

        self.factory.addPoint(R_LE[1, 0] * np.sin(-0.5 * theta), R_LE[1, 0] * np.cos(-0.5 * theta), x_inlet, 1.0, i_point)
        points_perio_2.append(i_point)
        i_point += 1

        for i in range(np.shape(X_LE)[1]):
            self.factory.addPoint(R_LE[1, i] * np.sin(0.5 * theta), R_LE[1, i] * np.cos(0.5 * theta), X_LE[1, i], 1.0,
                             i_point)
            points_perio_1.append(i_point)
            i_point += 1

            self.factory.addPoint(R_LE[1, i] * np.sin(-0.5 * theta), R_LE[1, i] * np.cos(-0.5 * theta), X_LE[1, i], 1.0,
                             i_point)
            points_perio_2.append(i_point)
            i_point += 1

            self.factory.addPoint(R_TE[1, i] * np.sin(0.5 * theta), R_TE[1, i] * np.cos(0.5 * theta), X_TE[1, i], 1.0,
                             i_point)
            points_perio_1.append(i_point)
            i_point += 1

            self.factory.addPoint(R_TE[1, i] * np.sin(-0.5 * theta), R_TE[1, i] * np.cos(-0.5 * theta), X_TE[1, i], 1.0,
                             i_point)
            points_perio_2.append(i_point)
            i_point += 1

        self.factory.addPoint(R_LE[1, 0] * np.sin(0.5 * theta), R_LE[1, 0] * np.cos(0.5 * theta), x_outlet, 1.0, i_point)
        points_perio_1.append(i_point)
        i_point += 1

        self.factory.addPoint(R_LE[1, 0] * np.sin(-0.5 * theta), R_LE[1, 0] * np.cos(-0.5 * theta), x_outlet, 1.0, i_point)
        points_perio_2.append(i_point)

        self.points_perio_1 = points_perio_1
        self.points_perio_2 = points_perio_2