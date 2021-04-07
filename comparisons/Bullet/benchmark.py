import json
import argparse
import pathlib
import datetime
import subprocess

import numpy
import tqdm
import pybullet as bullet
import igl


def get_time_stamp():
    return datetime.datetime.now().strftime("%Y-%b-%d-%H-%M-%S")


def save_mesh(meshes, out_path, index):
    Vs = []
    Fs = []
    offset = 0

    for i in range(bullet.getNumBodies()):
        try:
            V = meshes[i][0].copy()
        except:
            breakpoint()
        F = meshes[i][1].copy()

        t, q = bullet.getBasePositionAndOrientation(i)
        R = numpy.array(bullet.getMatrixFromQuaternion(q)).reshape(3, 3)

        V = V @ R.T + t

        F += offset
        offset += V.shape[0]

        Vs.append(V)
        Fs.append(F)

    V = numpy.concatenate(Vs)
    F = numpy.concatenate(Fs)

    igl.write_triangle_mesh(str(out_path / f"m_{index:04d}.obj"), V, F)


def object_from_obj(filename, mass=1, mesh_scale=[1, 1, 1]):
    col_flags = bullet.URDF_INITIALIZE_SAT_FEATURES
    if mass == 0:  # static objects
        col_flags = bullet.GEOM_FORCE_CONCAVE_TRIMESH

    collision_shape_id = bullet.createCollisionShape(
        shapeType=bullet.GEOM_MESH,
        fileName=str(filename),
        collisionFramePosition=[0, 0, 0],
        meshScale=mesh_scale,
        flags=col_flags)

    visual_shape_id = bullet.createVisualShape(
        shapeType=bullet.GEOM_MESH,
        fileName=str(filename),
        visualFramePosition=[0, 0, 0],
        meshScale=mesh_scale)

    body = bullet.createMultiBody(
        baseMass=mass,
        baseInertialFramePosition=[0, 0, 0],
        baseCollisionShapeIndex=collision_shape_id,
        baseVisualShapeIndex=visual_shape_id,
        basePosition=[0, 0, 0])

    average_vertex = numpy.zeros(3)
    if mass > 0:
        # compute inertial frame, assuming mass at vertices
        average_vertex = numpy.average(
            numpy.array(bullet.getMeshData(body)[1]), axis=0)

        # bullet.resetBasePositionAndOrientation(body, [100000,0,0],[0,0,0,1])
        bullet.removeBody(body)
        shift = -average_vertex

        collision_shape_id2 = bullet.createCollisionShape(
            shapeType=bullet.GEOM_MESH,
            fileName=str(filename),
            collisionFramePosition=shift,
            meshScale=mesh_scale,
            flags=col_flags)
        visual_shape_id2 = bullet.createVisualShape(
            shapeType=bullet.GEOM_MESH,
            fileName=str(filename),
            visualFramePosition=shift,
            meshScale=mesh_scale)

        body = bullet.createMultiBody(
            baseMass=mass,
            baseInertialFramePosition=[0, 0, 0],
            baseCollisionShapeIndex=collision_shape_id2,
            baseVisualShapeIndex=visual_shape_id2,
            basePosition=[0, 0, 0])

        # make dynamics objects random color
        color = (numpy.random.random(3)).tolist() + [1]

        bullet.changeVisualShape(body, -1, rgbaColor=color)

    bullet.changeVisualShape(
        body, -1, flags=bullet.VISUAL_SHAPE_DOUBLE_SIDED)

    return body, average_vertex


def convert_to_convex_mesh(in_path, out_path):
    bullet.vhacd(
        str(in_path), str(out_path), "vhacd_log.txt", concavity=0,
        maxNumVerticesPerCH=1024, depth=32, resolution=1000000,
        convexhullApproximation=0)


def run_simulation(fixture, meshes_path, out_path, args):
    rigid_body_problem = fixture["rigid_body_problem"]

    gravity = rigid_body_problem.get("gravity", [0, 0, 0])

    timestep = args.timestep
    if timestep is None:
        timestep = fixture.get("timestep", 1e-2)

    bullet.connect(bullet.DIRECT)
    bullet.setPhysicsEngineParameter(
        enableSAT=args.enable_sat,
        enableConeFriction=args.enable_cone_friction)
    bullet.setGravity(*gravity)
    bullet.setTimeStep(timestep)

    # plane = bullet.loadURDF("plane_implicit.urdf")#, pos = object_from_obj("meshes/plane.obj", mass=0)
    # bullet.changeDynamics(plane,-1,lateralFriction=1, frictionAnchor = 1, contactStiffness=30000, contactDamping=10000)
    # orn = bullet.getQuaternionFromEuler([ -numpy.pi/2,0, 0])
    # bullet.resetBasePositionAndOrientation(plane, [0,0,0], orn)

    meshes = []

    convex_meshes_path = pathlib.Path(__file__).resolve().parent / "meshes"

    for body in rigid_body_problem["rigid_bodies"]:
        if not body.get("enabled", True):
            continue

        mesh_path = meshes_path / body["mesh"]

        mass = body.get("density", 1000) # Assumes the volume is 1 m³

        is_dof_fixed = body.get("is_dof_fixed", False)
        if (body.get("type", "dynamic") == "static"
                or (isinstance(is_dof_fixed, list) and all(is_dof_fixed))
                or is_dof_fixed):
            mass = 0  # static object

        if args.make_convex and mass > 0:
            try:
                convex_mesh_path = (
                    convex_meshes_path / mesh_path.relative_to(meshes_path))
                convex_mesh_path.parent.mkdir(parents=True, exist_ok=True)
            except:
                convex_mesh_path = convex_meshes_path / mesh_path.name
            if not convex_mesh_path.exists():
                convert_to_convex_mesh(mesh_path, convex_mesh_path)
            mesh_path = convex_mesh_path

        V, F = igl.read_triangle_mesh(str(mesh_path))
        mesh_scale = body.get("scale", [1, 1, 1])
        if isinstance(mesh_scale, float):
            mesh_scale = [mesh_scale] * 3
        if "dimensions" in body:
            org_dim = V.max(axis=0) - V.min(axis=0)
            mesh_scale = body["dimensions"] / org_dim
        meshes.append((V * mesh_scale, F))

        body_id, com_shift = object_from_obj(
            mesh_path, mass=mass, mesh_scale=mesh_scale)

        # combined_friction = friction_a * friction_b so
        # keep friction of static objects to "default"=1
        mu = args.default_friction
        if mass > 0:
            mu = rigid_body_problem.get("coefficient_friction", 0.0)

        bullet.changeDynamics(
            body_id, -1, lateralFriction=mu, frictionAnchor=1)

        pos = body.get("position", [0, 0, 0])
        eul = numpy.deg2rad(body.get("rotation", [0, 0, 0]))
        orn = bullet.getQuaternionFromEuler(eul)

        lin_vel = body.get("linear_velocity", [0, 0, 0])
        ang_vel = numpy.deg2rad(body.get("angular_velocity", [0, 0, 0]))

        com_pos = pos + args.shift_mag * com_shift

        bullet.resetBasePositionAndOrientation(body_id, com_pos, orn)
        bullet.resetBaseVelocity(body_id, lin_vel, ang_vel)

    max_time = fixture.get("max_time", 5)

    index = 0
    prev_save = 0
    save_mesh(meshes, out_path, 0)

    for i in tqdm.tqdm(range(int(numpy.ceil(max_time / timestep)))):
        bullet.stepSimulation()
        if timestep >= 1e-2 or i * timestep - prev_save > 1e-2:
            index += 1
            save_mesh(meshes, out_path, index)
            prev_save = i * timestep


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test Rigid IPC examples in Bullet")
    parser.add_argument(
        "-i", "--input", "--json-file", metavar="path/to/input",
        type=pathlib.Path, dest="input", help="path to input json(s)",
        nargs="+")
    parser.add_argument(
        "--shift-mag", type=int, default=0,
        help=("Shift the collision/visual (due to OBJ file not centered) "
              "values -1=negative shift, 0=no shift, 1=positive shift"))
    parser.add_argument(
        "--camera-distance", type=float, default=7,
        help="Camera Distance (to target)")
    parser.add_argument(
        "--camera-target", type=float, default=[-1, 1, 1], nargs=3,
        help="Camera target (lookat) position")
    parser.add_argument(
        "--disable-sat", action="store_false", default=True, dest="enable_sat",
        help="Disable Separating Axis Test (SAT) collision")
    parser.add_argument(
        "--disable-cone-friction", action="store_false", default=True,
        dest="enable_cone_friction",
        help="Disable Cone friction (instead use the default pyramid friction model)")
    parser.add_argument(
        "--default-friction", type=float, default=0,
        help="Default friction coefficient (it nof provided in json file)")
    parser.add_argument(
        "--make-convex", action="store_true", default=False,
        help="Convert dynamic bodies to convex meshes (using V-HACD)")
    parser.add_argument(
        "--dt", "--timestep", type=float, default=None, dest="timestep",
        help="timestep")
    args = parser.parse_args()

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

    root_path = pathlib.Path(__file__).resolve().parents[2]
    meshes_path = root_path / "meshes"

    # renderer = root_path / "build" / "release" / "tools" / "render_simulation"
    # if not renderer.exists():
    renderer = None

    for input in args.input:
        try:
            out_path = input.relative_to(
                root_path / "fixtures" / "3D").with_suffix("")
        except:
            out_path = input.stem
        out_path = ("output" / out_path).resolve()
        if args.timestep is not None:
            out_path /= f"timestep={args.timestep:g}"
        print("out_path:", out_path)
        out_path.mkdir(exist_ok=True, parents=True)

        with open(input) as f:
            fixture = json.load(f)

        run_simulation(fixture, meshes_path, out_path, args)
        bullet.resetSimulation()

        # Render simulation
        if renderer is not None:
            print("Rendering simulation")
            video_name = f"{input.stem}-{get_time_stamp()}-chrono.mp4"
            subprocess.run([str(renderer),
                            "-i", out_path,
                            "-o", out_path / video_name,
                            "--fps", "100"])

        print()


if __name__ == "__main__":
    main()
