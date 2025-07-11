import copy
import numpy as np
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from ..utilities import print_dict

try:
    import gmsh

    GMSH_AVAILABLE = True
except ImportError:
    gmsh = None
    GMSH_AVAILABLE = False


def check_gmsh():
    if not GMSH_AVAILABLE:
        raise ImportError("gmsh is not installed. Run `pip install gmsh`.")


class NozzleMesher:
    def __init__(self, config):

        # Check gmsh is installed

        self.config = copy.deepcopy(config)

        # Validate upper side
        us = self.config["upper_side"]
        if us.get("type") != "explicit" or "segments" not in us:
            raise ValueError(
                "upper_side must be explicit and define one or several B-spline segments"
            )

        # Validate lower side
        ls = self.config["lower_side"]
        if ls.get("type") not in ("explicit", "mirror", "symmetry"):
            raise ValueError("lower_side.type must be explicit, mirror, or symmetry")
        if ls["type"] == "explicit" and "segments" not in ls:
            raise ValueError("lower_side explicit requires segments")

        # Unpack geometry and mesh settings
        self.lc = config["mesh"]["length_scale"]
        self.upper_segments = us["segments"]
        self.lower_type = ls.get("type", "symmetry")
        self.lower_segments = ls.get("segments", [])
        self.mesh_cfg = config["mesh"]
        self.Nx = self.mesh_cfg["Nx"]
        self.Ny = self.mesh_cfg["Ny"]
        self.has_inflation_layers = not np.isclose(
            self.mesh_cfg["inflation_ratio"], 1.0, atol=1e-6
        )

        # Storage for Gmsh tags
        self.curves = {}
        self.points = {}
        self.surfaces = {}
        self.corners = {}

        # Normalize lower segments (mirror, symmetry, or explicit)
        self.lower_segments = self._prepare_lower_segments()

        # Initialize Gmsh engine
        self._initialize()



    def create_mesh(self):
        """Create full nozzle mesh (geometry + structured meshing)."""
        
        # If not symmetry, then split the number of nodes among upper nad lower sides
        Ny_upper = int(np.ceil(self.Ny / 2)) if self.lower_type != "symmetry" else self.Ny
        Ny_lower = int(np.floor(self.Ny / 2)) if self.lower_type != "symmetry" else 0

        if not self.has_inflation_layers:
            
            # Create geometry
            self._build_upper_side(with_inflation=False)
            self._build_lower_side(with_inflation=False)
            
            # Upper core block
            self.surfaces["domain_upper"] = create_structured_surface(
                curves_1=self.curves["centerline"],
                curves_2=self.curves["outflow_core_upper"],
                curves_3=list(reversed(self.curves["wall_upper"])),  # CCW
                curves_4=self.curves["inflow_core_upper"],
                corner_12=self.points["outflow_centerline"],
                corner_23=self.points["outflow_wall_upper"],
                corner_34=self.points["inflow_wall_upper"],
                corner_41=self.points["inflow_centerline"],
                N_horizontal=self.Nx,
                N_vertical=Ny_upper,
            )

            # Lower core block
            if self.lower_type != "symmetry":
                self.surfaces["domain_lower"] = create_structured_surface(
                    curves_1=self.curves["wall_lower"],
                    curves_2=self.curves["outflow_core_lower"],
                    curves_3=self.curves["centerline"],
                    curves_4=self.curves["inflow_core_lower"],
                    corner_12=self.points["outflow_wall_lower"],
                    corner_23=self.points["outflow_centerline"],
                    corner_34=self.points["inflow_centerline"],
                    corner_41=self.points["inflow_wall_lower"],
                    N_horizontal=self.Nx,
                    N_vertical=Ny_lower,
                )

        else:

            self._build_upper_side(with_inflation=True, N_vertical=Ny_upper)
            self._build_lower_side(with_inflation=True, N_vertical=Ny_lower)

            # Upper core block
            self.surfaces["domain_upper_core"] = create_structured_surface(
                curves_1=self.curves["centerline"],
                curves_2=self.curves["outflow_core_upper"],
                curves_3=self.curves["inflation_upper"], 
                curves_4=self.curves["inflow_core_upper"],
                corner_12=self.points["outflow_centerline"],
                corner_23=self.points["outflow_inflation_upper"],
                corner_34=self.points["inflow_inflation_upper"],
                corner_41=self.points["inflow_centerline"],
                N_horizontal=self.Nx,
                N_vertical=self.inflation_data_upper["N_core"],
                progression_horizontal=1,
                progression_vertical=1,
            )

            # Upper inflation block
            self.surfaces["domain_upper_inflation"] = create_structured_surface(
                curves_1=self.curves["inflation_upper"],
                curves_2=self.curves["outflow_inflation_upper"],
                curves_3=list(reversed(self.curves["wall_upper"])),  # CCW
                curves_4=self.curves["inflow_inflation_upper"],
                corner_12=self.points["outflow_inflation_upper"],
                corner_23=self.points["outflow_wall_upper"],
                corner_34=self.points["inflow_wall_upper"],
                corner_41=self.points["inflow_inflation_upper"],
                N_horizontal=self.Nx,
                N_vertical=self.inflation_data_upper["N_inflation"],
                progression_horizontal=1.0,
                progression_vertical=-self.inflation_data_upper["inflation_ratio"],
            )

            if self.lower_type != "symmetry":
                # Lower core block
                self.surfaces["domain_lower_core"] = create_structured_surface(
                    curves_1=self.curves["centerline"],
                    curves_2=self.curves["outflow_core_lower"],
                    curves_3=self.curves["inflation_lower"], 
                    curves_4=self.curves["inflow_core_lower"],
                    corner_12=self.points["outflow_centerline"],
                    corner_23=self.points["outflow_inflation_lower"],
                    corner_34=self.points["inflow_inflation_lower"],
                    corner_41=self.points["inflow_centerline"],
                    N_horizontal=self.Nx,
                    N_vertical=self.inflation_data_lower["N_core"],
                    progression_horizontal=1,
                    progression_vertical=1,
                )

                # Lower inflation block
                self.surfaces["domain_lower_inflation"] = create_structured_surface(
                    curves_1=self.curves["inflation_lower"],
                    curves_2=self.curves["outflow_inflation_lower"],
                    curves_3=list(reversed(self.curves["wall_lower"])),  # CCW
                    curves_4=self.curves["inflow_inflation_lower"],
                    corner_12=self.points["outflow_inflation_lower"],
                    corner_23=self.points["outflow_wall_lower"],
                    corner_34=self.points["inflow_wall_lower"],
                    corner_41=self.points["inflow_inflation_lower"],
                    N_horizontal=self.Nx,
                    N_vertical=self.inflation_data_lower["N_inflation"],
                    progression_horizontal=1.0,
                    progression_vertical=-self.inflation_data_lower["inflation_ratio"],
                )

        # Create the structured mesh
        gmsh.model.mesh.generate(2)

        # # Extrude all fluid surfaces to create a 3D mesh
        # extruded_volumes = []
        # for name, surf_tag in self.surfaces.items():
        #     if "domain" in name:
        #         # Extrude 1.0 m in Z-direction, structured hex mesh
        #         out = gmsh.model.occ.extrude([(2, surf_tag)], 0, 0, 1.0, numElements=[1], recombine=True)
        #         extruded_volumes.extend(out)

        # gmsh.model.occ.synchronize()
        # gmsh.model.mesh.generate(3)


        self.create_physical_groups()
                  
    
    def write_mesh(self, filename="nozzle.msh"):
        # gmsh.option.setNumber("Mesh.MshFileVersion", 4.2)
        # gmsh.option.setNumber("Mesh.Binary", 0)  # ASCII output
        gmsh.write(filename)



    def show_mesh(self):
        black = [0, 0, 0]
        # Always check what is defined
        for name in [
            "domain_upper", "domain_lower",
            "domain_upper_core", "domain_lower_core",
            "domain_upper_inflation", "domain_lower_inflation"
        ]:
            if name in self.surfaces:
                gmsh.model.setColor([(2, self.surfaces[name])], *black)
                
        gmsh.fltk.run()
        # gmsh.finalize()

    def _initialize(self):
        """Initialize the Gmsh model."""
        check_gmsh()
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add(self.config["model_name"])


    def _build_upper_side(self, with_inflation: bool, N_vertical: int = None):
        """
        Build the upper half geometry (with optional inflation region).

        Parameters
        ----------
        with_inflation : bool
            Whether to include inflation layer geometry between wall and core.
        N_vertical : int, optional
            Number of cells in the vertical direction (required if with_inflation=True).
        """
        # 1. Upper wall
        self._build_streamwise_curves(self.upper_segments, "wall_upper")

        # 2. Centerline
        x_last = self.upper_segments[-1]["control_points"][-1][0]
        centerline = [{"order": 1, "control_points": [(0.0, 0.0), (x_last, 0.0)]}]
        self._build_streamwise_curves(centerline, "centerline")

        if not with_inflation:
            # Direct connection between center and upper wall
            self.curves["inflow_core_upper"] = [
                gmsh.model.occ.addLine(self.points["inflow_centerline"], self.points["inflow_wall_upper"])
            ]
            self.curves["outflow_core_upper"] = [
                gmsh.model.occ.addLine(self.points["outflow_centerline"], self.points["outflow_wall_upper"])
            ]
            gmsh.model.occ.synchronize()
            return

        # --- Inflation layer logic ---
        assert N_vertical is not None, "N_vertical must be provided when with_inflation=True"

        # 3. Generate 1D mesh for node extraction
        gmsh.model.occ.synchronize()
        set_transfinite_curves(self.curves["wall_upper"], self.Nx + 1)
        set_transfinite_curves(self.curves["centerline"], self.Nx + 1)
        gmsh.model.mesh.generate(1)

        # 4. Compute inflation distribution
        _, y, _ = gmsh.model.getValue(0, self.points["inflow_wall_upper"], [])
        self.inflation_data_upper = find_optimal_inflation_layer(
            H_total=y,
            inflation_ratio=self.mesh_cfg["inflation_ratio"],
            h_wall=self.mesh_cfg["first_cell_height"],
            N_total=N_vertical,
        )
        v = self.inflation_data_upper["H_core"] / self.inflation_data_upper["H_total"]

        # 5. Build B-spline curve for the inflation interface
        center_pts = sample_meshed_curves(self.curves["centerline"])
        upper_pts = sample_meshed_curves(self.curves["wall_upper"])
        CPs = [( (1 - v) * x0 + v * x1, (1 - v) * y0 + v * y1 )
            for (x0, y0, _), (x1, y1, _) in zip(center_pts, upper_pts)]
        segments = [{"order": 3, "control_points": CPs}]
        self._build_streamwise_curves(segments, "inflation_upper")

        # 6. Create vertical connections
        self.curves["inflow_core_upper"] = [
            gmsh.model.occ.addLine(self.points["inflow_centerline"], self.points["inflow_inflation_upper"])
        ]
        self.curves["outflow_core_upper"] = [
            gmsh.model.occ.addLine(self.points["outflow_centerline"], self.points["outflow_inflation_upper"])
        ]
        self.curves["inflow_inflation_upper"] = [
            gmsh.model.occ.addLine(self.points["inflow_inflation_upper"], self.points["inflow_wall_upper"])
        ]
        self.curves["outflow_inflation_upper"] = [
            gmsh.model.occ.addLine(self.points["outflow_inflation_upper"], self.points["outflow_wall_upper"])
        ]

        gmsh.model.occ.synchronize()

    def _build_lower_side(self, with_inflation: bool, N_vertical: int = None):
        """
        Build the lower half geometry (with optional inflation region).

        Parameters
        ----------
        with_inflation : bool
            Whether to include inflation layer geometry between wall and core.
        N_vertical : int, optional
            Number of cells in the vertical direction (required if with_inflation=True).
        """
        if self.lower_type == "symmetry":
            return  # No lower wall to define

        # 1. Lower wall
        self._build_streamwise_curves(self.lower_segments, "wall_lower")

        if not with_inflation:
            # Direct connection between center and lower wall
            self.curves["inflow_core_lower"] = [
                gmsh.model.occ.addLine(self.points["inflow_wall_lower"], self.points["inflow_centerline"])
            ]
            self.curves["outflow_core_lower"] = [
                gmsh.model.occ.addLine(self.points["outflow_wall_lower"], self.points["outflow_centerline"])
            ]
            gmsh.model.occ.synchronize()
            return

        # --- Inflation layer logic ---
        assert N_vertical is not None, "N_vertical must be provided when with_inflation=True"

        # 2. Generate 1D mesh for node extraction
        gmsh.model.occ.synchronize()
        set_transfinite_curves(self.curves["wall_lower"], self.Nx + 1)
        set_transfinite_curves(self.curves["centerline"], self.Nx + 1)
        gmsh.model.mesh.generate(1)

        # 3. Compute inflation distribution
        _, y, _ = gmsh.model.getValue(0, self.points["inflow_wall_lower"], [])
        self.inflation_data_lower = find_optimal_inflation_layer(
            H_total=abs(y),
            inflation_ratio=self.mesh_cfg["inflation_ratio"],
            h_wall=self.mesh_cfg["first_cell_height"],
            N_total=N_vertical,
        )
        v = self.inflation_data_lower["H_core"] / self.inflation_data_lower["H_total"]

        # 4. Build B-spline curve for the inflation interface (below centerline)
        center_pts = sample_meshed_curves(self.curves["centerline"])
        lower_pts = sample_meshed_curves(self.curves["wall_lower"])
        CPs = [( (1 - v) * x0 + v * x1, (1 - v) * y0 + v * y1 )
            for (x0, y0, _), (x1, y1, _) in zip(center_pts, lower_pts)]
        segments = [{"order": 3, "control_points": CPs}]
        self._build_streamwise_curves(segments, "inflation_lower")

        # 5. Create vertical connections
        self.curves["inflow_core_lower"] = [
            gmsh.model.occ.addLine(self.points["inflow_centerline"], self.points["inflow_inflation_lower"])
        ]
        self.curves["outflow_core_lower"] = [
            gmsh.model.occ.addLine(self.points["outflow_centerline"], self.points["outflow_inflation_lower"])
        ]
        self.curves["inflow_inflation_lower"] = [
            gmsh.model.occ.addLine(self.points["inflow_inflation_lower"], self.points["inflow_wall_lower"])
        ]
        self.curves["outflow_inflation_lower"] = [
            gmsh.model.occ.addLine(self.points["outflow_inflation_lower"], self.points["outflow_wall_lower"])
        ]

        gmsh.model.occ.synchronize()


    def _build_streamwise_curves(self, segments, name):
        """
        Build a list of streamwise curves from a sequence of B-spline segments.

        Each segment is a dictionary with:
            - 'order': degree of the B-spline (int)
            - 'control_points': list of (x, y) tuples

        If two consecutive segments share a control point at the junction
        (within a small absolute tolerance), the OCC point is reused
        to ensure geometric continuity and avoid duplication.

        Parameters
        ----------
        segments : list of dict
            Wall definition as a list of spline segments.

        name : str
            Base name under which to store the generated curve tags and point tags
            in self.curves and self.points, respectively.

        Returns
        -------
        (first_point_tag, last_point_tag) : tuple of int
            OCC point tags of the first and last control points in the wall.
        """

        def points_close(p1, p2, tol=1e-9):
            return abs(p1[0] - p2[0]) < tol and abs(p1[1] - p2[1]) < tol

        curve_tags = []
        prev_pt_tag = None
        prev_coords = None
        first_pt_tag = None

        for i, segment in enumerate(segments):
            ctrl_pts = segment["control_points"]
            degree = segment.get("order", len(ctrl_pts) - 1)

            gmsh_pts = []

            for j, (x, y) in enumerate(ctrl_pts):
                if i > 0 and j == 0 and points_close((x, y), prev_coords):
                    # Reuse point at segment junction if coordinates match
                    gmsh_pts.append(prev_pt_tag)
                else:
                    pt_tag = gmsh.model.occ.addPoint(x, y, 0.0, self.lc)
                    gmsh_pts.append(pt_tag)

            # Add B-spline curve for the current segment
            curve_tag = gmsh.model.occ.addBSpline(gmsh_pts, degree=degree)
            curve_tags.append(curve_tag)

            # Store first and last point tags
            if i == 0:
                first_pt_tag = gmsh_pts[0]
            prev_pt_tag = gmsh_pts[-1]
            prev_coords = ctrl_pts[-1]

        # Store results in geometry containers
        self.curves[name] = curve_tags
        self.points[f"inflow_{name}"] = first_pt_tag
        self.points[f"outflow_{name}"] = prev_pt_tag


    def _prepare_lower_segments(self):
        """
        Generate the list of B-spline segment definitions for the lower wall.

        Returns
        -------
        list of dict
            Each segment is a dictionary with keys:
                - 'order': spline order (int)
                - 'control_points': list of (x, y) tuples

        Behavior depends on self.lower_type:
        - "explicit" : Use self.lower_segments directly.
        - "mirror"   : Mirror upper segments across y=0 (negate y of control points).
        - "symmetry" : Generate a single straight segment along y=0 spanning x-range of upper wall.
        """
        if self.lower_type == "explicit":
            # Use explicitly defined lower wall segments
            return self.lower_segments

        if self.lower_type == "mirror":
            # Mirror upper wall across the x-axis (negate y)
            mirrored = []
            for seg in self.upper_segments:
                cps = [(x, -y) for x, y in seg["control_points"]]
                mirrored.append(
                    {"order": seg.get("order", len(cps) - 1), "control_points": cps}
                )
            return mirrored

        if self.lower_type == "symmetry":
            # Create a single flat lower wall from (0, 0) to (x_last, 0)
            x_last = self.upper_segments[-1]["control_points"][-1][0]
            return [{"order": 1, "control_points": [(0.0, 0.0), (x_last, 0.0)]}]

        raise ValueError(f"Unknown lower_type: {self.lower_type}")



    def create_physical_groups(self):
        """
        Define physical groups for Fluent boundary conditions based on self.curves.

        Parameters
        ----------
        symmetry : bool
            Whether the geometry includes a centerline symmetry boundary.
        """
        inlet_lines = []
        outlet_lines = []
        wall_lines = []
        symmetry_lines = []

        for name, tags in self.curves.items():
            lname = name.lower()
            if "inflow" in lname:
                inlet_lines.extend(tags)
            elif "outflow" in lname:
                outlet_lines.extend(tags)
            elif "wall" in lname:
                wall_lines.extend(tags)
            elif "centerline" in lname:
                symmetry_lines.extend(tags)

        if inlet_lines:
            pg = gmsh.model.addPhysicalGroup(1, inlet_lines)
            gmsh.model.setPhysicalName(1, pg, "inlet")

        if outlet_lines:
            pg = gmsh.model.addPhysicalGroup(1, outlet_lines)
            gmsh.model.setPhysicalName(1, pg, "outlet")

        if wall_lines:
            pg = gmsh.model.addPhysicalGroup(1, wall_lines)
            gmsh.model.setPhysicalName(1, pg, "wall")

        if self.lower_type == "symmetry":
            if not symmetry_lines:
                raise RuntimeError("Symmetry is enabled but no 'centerline' boundary found in self.curves.")
            pg = gmsh.model.addPhysicalGroup(1, symmetry_lines)
            gmsh.model.setPhysicalName(1, pg, "symmetry")

        # Handle fluid domain tagging
        fluid_surfaces = []
        for name, tags in self.surfaces.items():
            if "domain" in name.lower():
                fluid_surfaces.append(tags)

        if fluid_surfaces:
            pg = gmsh.model.addPhysicalGroup(2, fluid_surfaces)
            gmsh.model.setPhysicalName(2, pg, "fluid")
        else:
            raise RuntimeError("No surface with 'domain' found in self.surfaces. Cannot define 'fluid' region.")


def set_transfinite_curves(curve_tags, total_nodes, meshType="Progression", coef=1.0):
    """
    Assign transfinite meshing to a sequence of connected curves,
    distributing the given number of *points* (not elements) according to the
    arc length of each segment.

    Parameters:
    ----------
    curve_tags : list of int
        Tags of the OCC curves forming a continuous chain (e.g., lower or upper wall).
        Each curve is treated as a segment between two mesh nodes.

    total_nodes : int
        Total number of mesh nodes along the full chain (including endpoints).
        This implies total_nodes - 1 mesh elements will be distributed.

    Behavior:
    --------
    - If the chain has only one curve, assigns total_nodes directly to it.
    - If the chain has multiple curves, it:
        * Computes each curve’s arc length.
        * Distributes total_nodes - 1 elements proportionally to segment length.
        * Rounds to integers and ensures the sum matches exactly.
        * Applies `setTransfiniteCurve` to each segment with (elements + 1) points.
    """
    num_curves = len(curve_tags)

    if total_nodes < num_curves + 1:
        raise ValueError(
            f"total_nodes={total_nodes} must be at least number of segments+1 ({num_curves + 1})"
        )

    # Handle simple case: only one curve in the chain
    if num_curves == 1:
        gmsh.model.mesh.setTransfiniteCurve(
            curve_tags[0], total_nodes, meshType=meshType, coef=coef
        )
        return

    # Step 1: compute arc lengths for all segments
    lengths = np.array([get_arclength_occ(tag) for tag in curve_tags], dtype=float)
    total_length = lengths.sum()
    if total_length <= 0:
        raise ValueError("curve chain has zero total length")

    # Step 2: distribute elements proportionally to arc length
    total_elems = total_nodes - 1
    raw_elems = total_elems * lengths / total_length
    elems = np.maximum(1, np.round(raw_elems)).astype(int)
    # ensure at least 1 per segment

    # Step 3: correct any rounding error to exactly match total_elems
    elems[-1] += total_elems - elems.sum()
    assert elems.sum() == total_elems, "element count mismatch after correction"

    # Step 4: assign transfinite points for each curve (elements + 1 points)
    for tag, n_elem in zip(curve_tags, elems):
        gmsh.model.mesh.setTransfiniteCurve(
            tag, n_elem + 1, meshType=meshType, coef=coef
        )


def sample_meshed_curves(curve_tags):
    """
    Extracts ordered XYZ points from one or more curve tags,
    avoiding duplicated nodes at segment joints.

    Parameters
    ----------
    curve_tags : int or list of int
        A single curve tag or a list of curve tags.

    Returns
    -------
    list of tuple
        Ordered list of (x, y, z) points along the concatenated curves.
    """
    # ensure we have a list
    if isinstance(curve_tags, int):
        curve_tags = [curve_tags]
    pts = []
    for i, c in enumerate(curve_tags):
        # get mesh nodes for this 1D curve
        node_tags, coords, params = gmsh.model.mesh.getNodes(dim=1, tag=c)
        # sort by param to get ordering along curve
        triples = sorted(
            zip(node_tags, zip(coords[0::3], coords[1::3], coords[2::3]), params),
            key=lambda t: t[2],
        )
        local_pts = [xyz for _, xyz, _ in triples]

        # get start and end vertex coords
        bnds = gmsh.model.getBoundary([(1, c)], oriented=True)
        _, start_tag = bnds[0]
        _, end_tag = bnds[1]
        x0, y0, z0 = gmsh.model.getValue(0, start_tag, [])
        xf, yf, zf = gmsh.model.getValue(0, end_tag, [])

        if i == 0:
            # add the first curve’s start point
            pts.append((x0, y0, z0))
        # append interior points and curve end
        pts.extend(local_pts)
        pts.append((xf, yf, zf))

    return pts


def create_structured_surface(
    curves_1: list[int],
    curves_2: list[int],
    curves_3: list[int],
    curves_4: list[int],
    corner_12: int,
    corner_23: int,
    corner_34: int,
    corner_41: int,
    N_horizontal: int,
    N_vertical: int,
    progression_horizontal=1.0,
    progression_vertical=1.0,
) -> int:
    """
    Create and mesh a structured quadrilateral surface with transfinite constraints.

    Parameters
    ----------
    name : str
        Optional surface name for user reference or debugging.

    curves_1 to curves_4 : list of int
        OCC curve tags forming the four sides of the surface,
        ordered counter-clockwise: bottom, right, top, left.

    corner_12, corner_23, corner_34, corner_41 : int
        OCC point tags at the surface corners in CCW order.
        These are used for transfinite meshing.

    N_horizontal : int
        Number of mesh divisions along curves_1 and curves_3.

    N_vertical : int
        Number of mesh divisions along curves_2 and curves_4.

    Returns
    -------
    int
        The Gmsh surface tag of the created structured surface.
    """
    loop = curves_1 + curves_2 + curves_3 + curves_4
    loop_tag = gmsh.model.occ.addCurveLoop(loop)
    surf_tag = gmsh.model.occ.addPlaneSurface([loop_tag])
    gmsh.model.occ.synchronize()

    set_transfinite_curves(curves_1, N_horizontal + 1, coef=progression_horizontal)
    set_transfinite_curves(curves_3, N_horizontal + 1, coef=progression_horizontal)
    set_transfinite_curves(curves_2, N_vertical + 1, coef=progression_vertical)
    set_transfinite_curves(curves_4, N_vertical + 1, coef=progression_vertical)

    corners = [corner_12, corner_23, corner_34, corner_41]
    gmsh.model.mesh.setTransfiniteSurface(surf_tag, cornerTags=corners)
    gmsh.model.mesh.setRecombine(2, surf_tag)

    return surf_tag


def get_arclength_occ(curveTag):
    """
    Compute the arc length of a Gmsh OCC curve using numerical integration.

    Parameters
    ----------
    curveTag : int
        Tag of the OCC curve entity (dimension 1) whose arc length is to be computed.

    Returns
    -------
    length : float
        Total arc length of the curve, computed by integrating the norm of the
        first derivative (dx/dt, dy/dt, dz/dt) over the curve's parametric domain.

    Notes
    -----
    - This function uses `gmsh.model.getDerivative` to obtain the first derivative
      of the OCC curve with respect to its parameter `u`, and numerically integrates
      the magnitude of the derivative using `scipy.integrate.quad`.
    - Requires that `gmsh.model.occ.synchronize()` has been called before use.
    """
    dim = 1  # OCC curves are dimension 1

    def integrand(u):
        values = gmsh.model.getDerivative(dim, curveTag, np.atleast_1d(u))
        arr = np.array(values).reshape(-1, 3)
        dx, dy, dz = arr[0]
        return np.sqrt(dx * dx + dy * dy + dz * dz)

    u1, u2 = gmsh.model.getParametrizationBounds(dim, curveTag)
    length, _ = scipy.integrate.quad(integrand, u1, u2)
    return length


def find_optimal_inflation_layer(
    h_wall, inflation_ratio, H_total, N_total, check_errors=True
):
    """
    Determine the optimal inflation layer configuration given a geometric growth rate
    and total wall-normal mesh height and resolution.

    The function finds the number of inflation layer cells such that the last inflation cell
    matches the height of the uniform core cells. It also recomputes the exact inflation
    ratio that achieves this match after rounding the number of inflation cells.

    Parameters
    ----------
    h_wall : float
        Height of the first cell at the wall (e.g., to meet y⁺ constraint).
    
    inflation_ratio : float
        Initial approximation for the geometric growth ratio (r > 1). If r ≈ 1, a uniform grid is used.

    H_total : float
        Total height of the mesh in the wall-normal direction.

    N_total : int
        Total number of mesh cells in the wall-normal direction.

    check_errors : bool, default=True
        Whether to check and report mismatch between core and inflation cell heights.

    Returns
    -------
    dict
        Dictionary containing resolved inflation and core region parameters:
        - "H_total": total height of the domain
        - "N_total": total number of cells
        - "h_wall": height of the first inflation cell
        - "N_inflation": integer number of inflation cells (after rounding)
        - "H_inflation": height of the inflation region
        - "h_inflation": height of the last inflation cell
        - "N_core": number of core cells
        - "H_core": height of the core region
        - "h_core": uniform height of core cells
        - "inflation_ratio": recomputed exact geometric ratio used

    Method
    ------
    1. Model the inflation layer as a geometric series:
         h_i = h_wall * r^(i - 1) for i = 1...N_inflation

    2. Compute the total inflation height:
         H_inflation = h_wall * (1 - r^N) / (1 - r)

    3. Impose a matching constraint:
         h_last = h_wall * r^(N - 1) = h_core = (H_total - H_inflation) / (N_total - N)

    4. Solve this constraint as a root-finding problem for a floating-point N_inflation.

    5. Round N_inflation to the nearest integer, then recompute the exact inflation ratio
       that matches the target H_inflation height.

    Notes
    -----
    - If inflation_ratio ≈ 1, the method falls back to uniform spacing
      with no dedicated inflation layer (N_inflation = 0).
    - This ensures a smooth transition between inflation and core regions,
      and that the inflation region integrates cleanly into the mesh.
    """

    r = inflation_ratio
    if abs(inflation_ratio - 1.0) < 1e-12:
        # --- Step 0: No inflation layer: uniform spacing
        check_errors = False
        converged = True
        N_inflation = 0
        h_inflation_last = 0.0
        H_inflation = 0.0
        N_core = N_total
        H_core = H_total
        h_core = H_total / N_total

    else:

        # --- Step 1: Define mismatch between inflation and core layer thickness ---
        def objective(N_inflation_guess):
            """Root function: match h_last (inflation) to h_core (uniform core)."""
            H_infl = h_wall * (1 - r**N_inflation_guess) / (1 - r)
            h_last = h_wall * r**(N_inflation_guess - 1)
            N_core = N_total - N_inflation_guess
            H_core = H_total - H_infl
            h_core = H_core / N_core
            return h_core - h_last

        # --- Step 2: Solve for floating-point N_inflation ---
        result = scipy.optimize.root_scalar(objective, x0=N_total - 1)
        converged = result.converged
        N_inflation_float = result.root

        # --- Step 3: Compute provisional last inflation height from floating-point result ---
        H_inflation = h_wall * (1 - r**N_inflation_float) / (1 - r)
        h_inflation_last = h_wall * r**(N_inflation_float - 1)
        N_core = N_total - N_inflation_float
        H_core = H_total - H_inflation
        h_core = H_core / N_core

        # --- Step 4: Round to nearest integer and recompute exact ratio ---
        N_inflation = int(np.ceil(N_inflation_float))

        def residual(r_corrected):
            """Fix r to match the rounded N_inflation and original H_inflation."""
            # return h_wall * (1 - r_eff**N_inflation) / (1 - r_eff) - H_inflation
            H_infl = h_wall * (1 - r_corrected**N_inflation) / (1 - r_corrected)
            h_last = h_wall * r_corrected**(N_inflation - 1)
            N_core = N_total - N_inflation
            H_core = H_total - H_infl
            h_core = H_core / N_core
            return h_core - h_last
        
        r = scipy.optimize.root_scalar(residual, x0=inflation_ratio).root

        # --- Step 5: Final geometry based on corrected r ---
        H_inflation = h_wall * (1 - r**N_inflation) / (1 - r)
        N_core = N_total - N_inflation
        H_core = H_total - H_inflation
        h_inflation_last = h_wall * r**(N_inflation - 1)
        h_core = H_core / N_core

    inflation_data = {
        # "r": r,  # skipped as requested
        "H_total": H_total,
        "N_total": N_total,
        "h_wall": h_wall,
        "N_inflation": N_inflation,
        "H_inflation": H_inflation,
        "h_inflation": h_inflation_last,
        "N_core": N_core,
        "h_core": h_core,
        "H_core": H_core,
        "inflation_ratio": r,
    }

    # Sanity checks on the computed results
    errors = []

    if not converged:
        errors.append(
            "Inflation layer calculation failed: the Newton method did not converge.\n"
        )

    if H_inflation > H_total:
        errors.append("H_inflation > H_total. Try reducing the inflation ratio.\n")

    if N_inflation < 1:
        errors.append(
            "The Newton method did not find a solution with a positive number of inflation cells.\n"
            "Try increasing the inflation ratio or the number of cells.\n"
        )

    relative_diff = abs(h_inflation_last - h_core) / h_core
    threshold = inflation_ratio - 1.0
    if relative_diff > threshold:
        errors.append(
            f"h_inflation and h_core differ by {relative_diff:.2%}, which exceeds the allowed threshold (r - 1 = {threshold:.2%}).\n"
            "This indicates a mismatch between the last inflation cell and the core mesh height.\n"
        )

    # Raise error with full report if any checks fail
    if errors and check_errors:
        msg = ["\n".join(errors), "Computed inflation layer parameters:"]
        for k, v in inflation_data.items():
            msg.append(
                f"\t{k:16}: {v:.4g}" if isinstance(v, float) else f"\t{k:16}: {v}"
            )
        raise ValueError("\n".join(msg))

    return inflation_data


def plot_y_mesh(h_wall, r, H_total, N_total):

    params = find_optimal_inflation_layer(h_wall, r, H_total, N_total)
    N_inflation = params["N_inflation"]
    N_core = params["N_core"]
    h_wall = params["h_wall"]
    r = params["inflation_ratio"]
    H = params["H_total"]

    # y positions of cell edges
    y_bl = [0.0]
    for i in range(N_inflation):
        h = h_wall * r**i
        y_bl.append(y_bl[-1] + h)
    y_bl = np.array(y_bl)

    y_core = np.linspace(y_bl[-1], H, N_core + 1)

    # Plot (vertical orientation)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(np.zeros_like(y_bl), y_bl, "ro-", label="Inflation layer")
    ax.plot(np.zeros_like(y_core), y_core, "bo-", label="Core region")

    for y in y_bl:
        ax.axhline(y, color="red", linestyle="-")
    for y in y_core:
        ax.axhline(y, color="blue", linestyle="--")

    # ax.set_ylim(0, H)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xticks([])
    ax.set_ylabel("y-position")
    ax.set_title(f"Y-direction mesh: {N_inflation} BL cells, {N_core} core cells")
    ax.legend(loc="upper right")
    plt.grid(False)
    plt.tight_layout(pad=1)
    plt.show()
