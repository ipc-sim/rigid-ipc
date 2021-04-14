# Script to test Rigid IPC examples in Project Chrono
# Adapted from a pychrono example originally from Alessandro Tasora

import json
import pathlib
import argparse
import subprocess
from datetime import datetime

import numpy
from scipy.spatial.transform import Rotation

import tqdm

import pychrono.core as chrono
import igl


timestepper_names = {
    0: "EULER_IMPLICIT_LINEARIZED",
    1: "EULER_IMPLICIT_PROJECTED",
    2: "EULER_IMPLICIT",
    3: "TRAPEZOIDAL",
    4: "TRAPEZOIDAL_LINEARIZED",
    5: "HHT",
    6: "HEUN",
    7: "RUNGEKUTTA45",
    8: "EULER_EXPLICIT",
    9: "LEAPFROG",
    10: "NEWMARK",
    20: "CUSTOM",
}

solver_names = {
    0: "PSOR",
    1: "PSSOR",
    2: "PJACOBI",
    3: "PMINRES",
    4: "BARZILAIBORWEIN",
    5: "APGD",
    6: "ADDM",
    7: "SPARSE_LU",
    8: "SPARSE_QR",
    9: "PARDISO_MKL",
    10: "PARDISO_PROJECT",
    11: "MUMPS",
    12: "GMRES",
    13: "MINRES",
    14: "BICGSTAB",
    15: "CUSTOM"
}

box_vertices = numpy.array([
    [-0.5, -0.5, -0.5],
    [0.5, 0.5, -0.5],
    [0.5, -0.5, -0.5],
    [-0.5, 0.5, -0.5],
    [-0.5, 0.5, 0.5],
    [-0.5, -0.5, 0.5],
    [0.5, 0.5, 0.5],
    [0.5, -0.5, 0.5]])
box_faces = numpy.array([
    [0, 1, 2],
    [0, 3, 1],
    [0, 4, 3],
    [0, 5, 4],
    [3, 6, 1],
    [3, 4, 6],
    [2, 1, 6],
    [2, 6, 7],
    [0, 2, 7],
    [0, 7, 5],
    [5, 7, 6],
    [5, 6, 4]], dtype=int)


def print_simulation_parameters(system):
    print(f"""SolverType: {solver_names[system.GetSolverType()]},
TimestepperType: {timestepper_names[system.GetTimestepperType()]},
SolverMaxIterations: {system.GetSolverMaxIterations():d},
SolverTolerance: {system.GetSolverTolerance():g},
Maxiter: {system.GetMaxiter():d},
MaxPenetrationRecoverySpeed: {system.GetMaxPenetrationRecoverySpeed():g},
DefaultSuggestedEnvelope: {chrono.ChCollisionModel.GetDefaultSuggestedEnvelope():g},
DefaultSuggestedMargin: {chrono.ChCollisionModel.GetDefaultSuggestedMargin():g},
""")


def set_simulation_parameters(system, collision_envelope=None, collision_margin=None):
    print("Default parameters:")
    print_simulation_parameters(system)

    # system.SetSolverType(chrono.ChSolver.Type_APGD)
    system.SetSolverType(chrono.ChSolver.Type_BARZILAIBORWEIN)
    system.SetTimestepperType(
        chrono.ChTimestepper.Type_EULER_IMPLICIT_PROJECTED)
    system.SetSolverMaxIterations(15000)
    system.SetSolverTolerance(1e-12)
    system.SetMaxiter(1000)
    # system.SetMaxPenetrationRecoverySpeed(0.1)

    # Set the global collision margins. This is expecially important for very large or
    # very small objects. Set this before creating shapes. Not before creating system.
    if collision_envelope is not None:
        chrono.ChCollisionModel.SetDefaultSuggestedEnvelope(collision_envelope)
    if collision_margin is not None:
        chrono.ChCollisionModel.SetDefaultSuggestedMargin(collision_margin)

    print("Using parameters:")
    print_simulation_parameters(system)


def get_time_stamp():
    return datetime.now().strftime("%Y-%b-%d-%H-%M-%S")


def get_trafos(body):
    t = body.GetPos()
    t = numpy.array([t.x, t.y, t.z])
    R = numpy.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            R[i, j] = body.GetA()[i, j]
    return t, R


def save_mesh(system, meshes, out_path, index):
    bds = system.Get_bodylist()

    Vs = []
    Fs = []
    offset = 0

    for bi, body in enumerate(bds):
        t, R = get_trafos(body)
        V = meshes[bi][0].copy()
        F = meshes[bi][1].copy()

        V = V @ R.T
        V += t

        F += offset
        offset += V.shape[0]

        Vs.append(V)
        Fs.append(F)

    V = numpy.concatenate(Vs)
    F = numpy.concatenate(Fs)

    igl.write_triangle_mesh(str(out_path / f"m_{index:04d}.obj"), V, F)


def run_simulation(fixture, mesh_path, out_path, timestep=None,
                   collision_envelope=None, collision_margin=None, mu=None):
    rb_problem = fixture["rigid_body_problem"]

    system = chrono.ChSystemNSC()
    system.Set_G_acc(chrono.ChVectorD(*rb_problem.get("gravity", [0, 0, 0])))
    dhat = fixture.get(
        "distance_barrier_constraint",
        {"initial_barrier_activation_distance": 1e-3}).get(
        "initial_barrier_activation_distance", 1e-3)
    if collision_envelope is not None and collision_envelope < 0:
        collision_envelope = dhat
    if collision_margin is not None and collision_margin < 0:
        collision_margin = dhat
    set_simulation_parameters(system, collision_envelope, collision_margin)

    # ---------------------------------------------------------------------
    #
    #  Create the simulation system and add items
    #

    # Create a contact material (with default properties, shared by all collision shapes)
    contact_material = chrono.ChMaterialSurfaceNSC()

    if mu is None:
        mu = rb_problem.get("coefficient_friction", 0)
    contact_material.SetFriction(max(mu, 0))
    # comp = system.GetMaterialCompositionStrategy()

    meshes = []

    for body in rb_problem["rigid_bodies"]:
        mpath = str(mesh_path / body["mesh"])

        if body["mesh"] == "plane.obj":
            pos = body["position"]
            dim = body.get("dimensions", [10, 0.005, 10])
            rot = body.get("rotation", [0, 0, 0])
            rot = Rotation.from_euler('xyz', rot, degrees=True)
            rot = rot.as_quat()

            floor = chrono.ChBodyEasyBox(*dim, 1, True, True, contact_material)
            floor.SetPos(chrono.ChVectorD(*pos))
            floor.SetRot(chrono.ChQuaternionD(rot[3], rot[0], rot[1], rot[2]))
            if body.get("is_dof_fixed", False) or body.get("type", None) == "static":
                floor.SetBodyFixed(True)
            system.Add(floor)
            meshes.append((dim * box_vertices, box_faces))
            continue

        if body.get("scale", 1) != 1:
            v, f = igl.read_triangle_mesh(mpath)
            v *= body["scale"]
            igl.write_triangle_mesh("tmp.obj", v, f)
            mpath = "tmp.obj"

        density = body.get("density", 1000)

        cbody = chrono.ChBodyEasyMesh(
            mpath,             # mesh filename
            density,           # density kg/m^3
            True,              # automatically compute mass and inertia?
            True,              # visualize?
            True,              # collide?
            contact_material)  # contact material

        cbody.Update()
        V, F = igl.read_triangle_mesh(mpath)
        t, R = get_trafos(cbody)
        # print(t)
        # print(R)
        V = V - t
        V = V @ R
        meshes.append((V, F))

        # rot = Rotation.from_euler('zyx', body["rotation"], degrees=True)
        rot2 = Rotation.from_matrix(R)
        rotTmp = body.get("rotation", [0.0, 0.0, 0.0])
        rot = Rotation.from_euler('xyz', rotTmp, degrees=True) * rot2
        rot = rot.as_quat()
        cbody.SetRot(chrono.ChQuaternionD(rot[3], rot[0], rot[1], rot[2]))

        if "position" in body:
            pos = body["position"]
            pos = [x + y for (x, y) in zip(pos, t)]
            cbody.SetPos(chrono.ChVectorD(pos[0], pos[1], pos[2]))

        # Need to fix coordinate systems
        if body.get("linear_velocity", [0, 0, 0]) != [0, 0, 0]:
            raise NotImplementedError(
                "Initial linear velocity not implemented!")
            # vel = body["linear_velocity"]
            # cbody.SetPos_dt(chrono.ChVectorD(vel[0], vel[1], vel[2]))
        if body.get("angular_velocity", [0, 0, 0]) != [0, 0, 0]:
            raise NotImplementedError(
                "Initial angular velocity not implemented!")
            # vel = body["angular_velocity"]
            # cbody.SetWVel_loc(chrono.ChVectorD(vel[0], vel[1], vel[2]))

        is_fixed = body.get("type", "dynamic") == "static"
        if not is_fixed:
            is_dof_fixed = body.get("is_dof_fixed", False)
            if isinstance(is_dof_fixed, list):
                is_fixed = all(is_dof_fixed)
                if not is_fixed and any(is_dof_fixed):
                    raise NotImplementedError("Only all dofs fixed supported!")
            else:
                is_fixed = is_dof_fixed
        cbody.SetBodyFixed(is_fixed)

        cbody.Update()
        system.Add(cbody)

    if timestep is None:
        timestep = fixture.get("timestep", 1e-2)
    max_time = fixture.get("max_time", 10)

    index = 0
    prev_save = 0
    save_mesh(system, meshes, out_path, 0)

    for i in tqdm.tqdm(range(int(numpy.ceil(max_time / timestep)))):
        system.DoStepDynamics(timestep)
        if timestep >= 1e-2 or system.GetChTime() - prev_save > 1e-2:
            index += 1
            save_mesh(system, meshes, out_path, index)
            prev_save = system.GetChTime()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test Rigid IPC examples in Project Chrono")
    parser.add_argument(
        "-i", "--input", metavar="path/to/input", type=pathlib.Path,
        dest="input", help="path to input json(s)", nargs="+")
    parser.add_argument(
        "--chrono-data", metavar="path/to/chrono/data", type=str,
        default="/usr/local/share/chrono/data/", dest="chrono_data_path",
        help="path to Chrono data")
    parser.add_argument(
        "--dt", "--timestep", type=float, default=None, dest="timestep",
        help="timestep")
    parser.add_argument(
        "--mu", type=float, default=None, dest="mu", help="coeff. friction")

    def env_mar_to_float(x):
        if x == "dhat":
            return -1
        return float(x)

    parser.add_argument("--envelope", type=env_mar_to_float, default=None,
                        help="collision envelope value")
    parser.add_argument("--margin", type=env_mar_to_float, default=None,
                        help="collision margin value")
    parser.add_argument("--no-video", action="store_true", default=False,
                        help="skip rendering")
    args = parser.parse_args()

    if not pathlib.Path(args.chrono_data_path).exists():
        parser.error(f"Invalid Chrono data path: {args.chrono_data_path}")

    inputs = []
    for input in args.input:
        if input.is_file() and input.suffix == ".json":
            inputs.append(input.resolve())
        elif input.is_dir():
            for glob_input in input.glob('**/*.json'):
                inputs.append(glob_input.resolve())
    args.input = inputs

    return args


def main():
    args = parse_args()

    # The path to the Chrono data directory containing various assets (meshes,
    # textures, data files) is automatically set, relative to the default
    # location of this demo. If running from a different directory, you must
    # change the path to the data directory with:
    chrono.SetChronoDataPath(str(args.chrono_data_path))

    root_path = pathlib.Path(__file__).resolve().parents[2]
    mesh_path = root_path / "meshes"

    renderer = root_path / "build" / "release" / "tools" / "render_simulation"
    if not renderer.exists():
        renderer = None

    cwd_output = pathlib.Path("output").resolve()

    for input in args.input:
        with open(input) as f:
            fixture = json.load(f)

        if args.timestep is None:
            args.timestep = fixture.get("timestep", 1e-2)

        try:
            out_path = input.resolve().relative_to(
                root_path / "fixtures" / "3D").with_suffix("")
        except:
            out_path = input.stem
        out_path = cwd_output / out_path
        folder_name = "_".join(
            ([] if args.timestep is None else [f"timestep={args.timestep:g}"])
            + ([] if args.mu is None else [f"mu={args.mu:g}"]))
        out_path /= folder_name
        print("out_path:", out_path)
        out_path.mkdir(exist_ok=True, parents=True)

        try:
            run_simulation(fixture, mesh_path, out_path,
                           timestep=args.timestep,
                           collision_envelope=args.envelope,
                           collision_margin=args.margin,
                           mu=args.mu)

            # Render simulation
            if renderer is not None and not args.no_video:
                print("Rendering simulation")
                video_name = f"{input.stem}-{get_time_stamp()}-chrono.mp4"
                subprocess.run([str(renderer),
                                "-i", out_path,
                                "-o", out_path / video_name,
                                "--fps", "100"])
        except NotImplementedError as err:
            print(f"{input}: {err}")
            out_path.rmdir()
        print()


if __name__ == "__main__":
    main()
