import meshio
import sys
import os

def convert_msh_to_dat(msh_file):
    """
    Converts a Gmsh .msh file to .dat format for Charlie.

    Outputs:
    - <base_name>.geo.dat  : Contains nodal coordinates and elements.
    - <base_name>.dom.dat  : Defines mesh dimensions and references geometry.
    - <base_name>.fix.dat  : Specifies boundary conditions.

    Parameters:
    msh_file (str): Path to the Gmsh-generated .msh file.

    Assumptions:
    - The mesh contains either TRI03 (2D) or TET04 (3D) elements.
    - Zero-based node indexing from Gmsh is corrected for Fortran (1-based).
    """

    # Ensure the file exists
    if not os.path.isfile(msh_file):
        print(f"❌ Error: Mesh file '{msh_file}' not found.")
        sys.exit(1)

    # Determine output directory and filenames
    output_dir = os.path.dirname(msh_file) or "."
    base_name = os.path.splitext(os.path.basename(msh_file))[0]  # Remove .msh extension
    geo_file = os.path.join(output_dir, f"{base_name}.geo.dat")
    dom_file = os.path.join(output_dir, f"{base_name}.dom.dat")
    fix_file = os.path.join(output_dir, f"{base_name}.fix.dat")

    print(f"📂 Processing mesh: {msh_file}")

    # Read the mesh
    mesh = meshio.read(msh_file)

    # Determine mesh properties
    num_nodes = len(mesh.points)
    num_elements = sum(len(c.data) for c in mesh.cells)
    num_boundaries = len(mesh.cells_dict.get("line", [])) + len(mesh.cells_dict.get("triangle", []))
    ndim = 3 if len(mesh.points[0]) == 3 else 2
    element_type = "TRI03" if "triangle" in mesh.cells_dict else "TET04"

    print(f"📄 Writing nodal coordinates to: {geo_file}")
    print(f"📄 Writing domain metadata to: {dom_file}")
    print(f"📄 Writing boundary conditions to: {fix_file}")

    # ✅ Generate `dom.dat` (Mesh Metadata)
    with open(dom_file, "w") as f:
        f.write("$--------------------------------------------------\n")
        f.write(" DIMENSIONS\n")
        f.write(f"   NODAL_POINTS    {num_nodes}\n")
        f.write(f"   ELEMENTS        {num_elements}\n")
        f.write(f"   BOUNDARIES      {num_boundaries}\n")
        f.write(f"   SPACE_DIMENSIONS= {ndim}\n")
        f.write(f"   TYPES_OF_ELEMS= {element_type}\n")
        f.write(" END_DIMENSIONS\n")
        f.write("$--------------------------------------------------\n")
        f.write(" GEOMETRY\n")
        f.write(f"   INCLUDE  {base_name}.geo.dat\n")
        f.write(" END_GEOMETRY\n")
        f.write("$--------------------------------------------------\n")
        f.write(" BOUNDARY_CONDITIONS\n")
        f.write("  ON_BOUNDARIES\n")
        f.write(f"   INCLUDE  {base_name}.fix.dat\n")
        f.write("  END_ON_BOUNDARIES\n")
        f.write("\n")
        f.write("  ON_NODES\n")
        f.write("  END_ON_NODES\n")
        f.write(" END_BOUNDARY_CONDITIONS\n")
        f.write("$--------------------------------------------------\n")

    # ✅ Generate `geo.dat` (Nodal Coordinates and Elements)
    with open(geo_file, "w") as f:
        f.write("COORDINATES\n")
        for i, node in enumerate(mesh.points, start=1):  # Convert to 1-based indexing
            f.write(f"{i} {node[0]} {node[1]}")
            if ndim == 3:
                f.write(f" {node[2]}")
            f.write("\n")
        f.write("END_COORDINATES\n")

        # 🔥 Add ELEMENTS section for connectivity (essential for 3D mesh!)
        if "triangle" in mesh.cells_dict or "tetra" in mesh.cells_dict:
            f.write("ELEMENTS\n")
            cell_type = "triangle" if ndim == 2 else "tetra"
            for i, cell in enumerate(mesh.cells_dict[cell_type], start=1):
                f.write(f"{i} " + " ".join(str(n + 1) for n in cell) + "\n")  # Convert to 1-based indexing
            f.write("END_ELEMENTS\n")

    # ✅ Generate `fix.dat` (Boundary Conditions)
    with open(fix_file, "w") as f:
        f.write("BOUNDARIES\n")
        boundary_found = False

        if "line" in mesh.cells_dict:  # 2D boundary edges
            for i, cell in enumerate(mesh.cells_dict["line"], start=1):
                f.write(f"{i} {cell[0]+1} {cell[1]+1} 1\n")  # Convert to 1-based indexing
                boundary_found = True

        if "triangle" in mesh.cells_dict:  # 3D boundary faces
            for i, cell in enumerate(mesh.cells_dict["triangle"], start=1):
                f.write(f"{i} {cell[0]+1} {cell[1]+1} {cell[2]+1} 1\n")  # Convert to 1-based indexing
                boundary_found = True

        if not boundary_found:
            print("⚠ Warning: No boundary conditions detected in the mesh.")

        f.write("END_BOUNDARIES\n")

    print("✅ Mesh conversion completed successfully!")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python utils/gmsh2charlie.py <mesh_file.msh>")
        sys.exit(1)

    msh_file = sys.argv[1]  # Get filename from command-line argument
    convert_msh_to_dat(msh_file)
