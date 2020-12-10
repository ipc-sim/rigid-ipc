import sys
import os
import pathlib
import json
import textwrap

import numpy
from scipy.spatial.transform import Rotation

import pymesh

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / "tools" / "fixtures"))  # noqa
import fixture_utils


def convert_to_ipc_msh(input_path, output_path):
    input_mesh = pymesh.load_mesh(str(input_path))
    if input_mesh.num_voxels == 0:  # tetrahedralize the mesh
        tetgen = pymesh.tetgen()
        tetgen.points = input_mesh.vertices
        tetgen.triangles = input_mesh.faces
        tetgen.split_boundary = False
        tetgen.max_radius_edge_ratio = 2.0
        tetgen.min_dihedral_angle = 0.0
        tetgen.merge_coplanar = False
        tetgen.run()
        output_mesh = tetgen.mesh
    else:  # mesh is already a tet mesh, just convert to IPC compatible
        output_mesh = input_mesh

    # pymesh.save_mesh(str(args.output), output_mesh, ascii=True)
    with open(output_path, mode='w') as f:
        f.write("$MeshFormat\n4 0 8\n$EndMeshFormat\n")
        f.write("$Entities\n0 0 0 1\n")
        f.write("0 {:g} {:g} {:g} {:g} {:g} {:g} 0 0\n".format(
            *output_mesh.nodes.min(axis=0), *output_mesh.nodes.max(axis=0)))
        f.write("$EndEntities\n")
        f.write("$Nodes\n")
        f.write("1 {0:d}\n0 3 0 {0:d}\n".format(output_mesh.num_nodes))
        for i, node in enumerate(output_mesh.nodes):
            f.write("{:d} {:g} {:g} {:g}\n".format(i + 1, *node))
        f.write("$EndNodes\n")

        f.write("$Elements\n")
        f.write("1 {0:d}\n0 3 4 {0:d}\n".format(output_mesh.num_elements))
        for i, element in enumerate(output_mesh.elements):
            f.write("{:d} {:d} {:d} {:d} {:d}\n".format(i + 1, *(element + 1)))
        f.write("$EndElements\n")

        f.write("$Surface\n")
        f.write("{:d}\n".format(output_mesh.num_faces))
        for face in output_mesh.faces:
            f.write("{0:d} {2:d} {1:d}\n".format(*(face + 1)))
        f.write("$EndSurface\n")


def fixture_to_ipc_script(fixture, output_path):
    timestep = fixture.get("timestep", 0.01)

    if "max_time" in fixture:
        max_time = fixture["max_time"]
    elif "max_iterations" in fixture:
        max_time = fixture["max_iterations"] * timestep
    else:
        max_time = 10

    dhat = fixture.get("distance_barrier_constraint", {}).get(
        "initial_barrier_activation_distance", 1e-3)

    bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    shapes = []
    for body in bodies:
        if "type" in body:
            is_static = body["type"] == "static"
        elif "is_dof_fixed" in body:
            is_static = (body["is_dof_fixed"] if isinstance(
                body["is_dof_fixed"], bool) else all(body["is_dof_fixed"]))
        else:
            is_static = False

        surface_mesh_path = pathlib.Path(body["mesh"])
        if not surface_mesh_path.exists():
            surface_mesh_path = (
                fixture_utils.get_meshes_dir_path() / surface_mesh_path)
        if not is_static:
            try:
                msh_path = surface_mesh_path.with_suffix('.msh').resolve()
                msh_path = (
                    pathlib.Path(__file__).resolve().parent / "meshes" /
                    msh_path.relative_to(fixture_utils.get_meshes_dir_path()))
                msh_path.parent.mkdir(parents=True, exist_ok=True)
            except:
                msh_path = surface_mesh_path.with_suffix('.msh')
            if not msh_path.exists():
                convert_to_ipc_msh(surface_mesh_path, msh_path)
            mesh_path = msh_path
        else:
            mesh_path = surface_mesh_path

        mesh_path = os.path.relpath(
            mesh_path.resolve(), output_path.parent.resolve())

        if "dimensions" in body:
            vertices = pymesh.load_mesh(str(surface_mesh_path)).vertices
            initial_dimensions = abs(
                vertices.max(axis=0) - vertices.min(axis=0))
            initial_dimensions[initial_dimensions == 0] = 1
            scale = numpy.array(body["dimensions"], dtype=float)
            scale /= initial_dimensions
        else:
            scale = body.get("scale", [1, 1, 1])
            if isinstance(scale, int) or isinstance(scale, float):
                scale = [scale] * 3

        if "rotation" in body:
            rotation = Rotation.from_euler(
                "xyz", body["rotation"], degrees=True).as_euler("zyx", degrees=True)[::-1]
        else:
            rotation = [0, 0, 0]

        linear_velocity = body.get("linear_velocity", [0, 0, 0])
        angular_velocity = body.get("angular_velocity", [0, 0, 0])
        force = body.get("force", [0, 0, 0])
        if "force" in body:
            nbc = ("  NBC -1e300 -1e300 -1e300  1e300 1e300 1e300 "
                   " {:g} {:g} {:g}".format(*body["force"]))
        else:
            nbc = ""
        if "torque" in body:
            print("External torque is not supported in IPC! Dropping it!")

        shapes.append(
            "{}  {:g} {:g} {:g}  {:g} {:g} {:g}  {:g} {:g} {:g} material {:g} 2e11 0.3  {}{}".format(
                mesh_path, *body.get("position", [0, 0, 0]),
                *rotation, *scale, body.get("density", 1000),
                "linearVelocity 0 0 0" if is_static else
                "initVel {:g} {:g} {:g}  {:g} {:g} {:g}".format(
                    *linear_velocity, *angular_velocity), nbc
            ))

    epsv = fixture.get("friction_constraints", {}).get(
        "static_friction_speed_bound", 1e-3)
    friction_iterations = fixture.get(
        "friction_constraints", {}).get("friction_iterations", 1)

    return textwrap.dedent(f"""\
        energy NH
        warmStart 0
        time {max_time} {timestep}

        shapes input {len(shapes)}
        {{}}

        selfCollisionOn
        selfFric {fixture["rigid_body_problem"].get("coefficient_friction", 0)}

        constraintSolver interiorPoint
        dHat {dhat}
        epsv {epsv}
        fricIterAmt {friction_iterations}
        """).format("\n".join(shapes))


def main():
    assert(len(sys.argv) > 1)
    input_path = pathlib.Path(sys.argv[1])

    try:
        output_path = input_path.with_suffix(".txt").resolve()
        output_path = (
            pathlib.Path(__file__).resolve().parent / "scripts" /
            output_path.relative_to(fixture_utils.get_fixture_dir_path()))
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except:
        output_path = input_path.with_suffix(".txt")

    with open(input_path) as input_file:
        fixture = json.load(input_file)

    ipc_script = fixture_to_ipc_script(fixture, output_path)

    with open(output_path, 'w') as ipc_script_file:
        ipc_script_file.write(ipc_script)


if __name__ == "__main__":
    main()
