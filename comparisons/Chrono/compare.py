# Script to test Rigid IPC examples in Project Chrono
# Adapted from a pychrono example originally from Alessandro Tasora


import pychrono.core as chrono
import pychrono.irrlicht as chronoirr
import time
import json
import os
import igl
import numpy as np
from scipy.spatial.transform import Rotation


def get_trafos(body):
    t = body.GetPos()
    t = np.array([t.x, t.y, t.z])
    R = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            R[i, j] = body.GetA()[i, j]

    return t, R


def save_mesh(application, out_path, index):
    sys = application.GetSystem()
    bds = sys.Get_bodylist()

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

    V = np.concatenate(Vs)
    F = np.concatenate(Fs)

    igl.write_triangle_mesh(os.path.join(
        out_path, "m_{:04d}.obj".format(index)), V, F)


if __name__ == "__main__":
    root_path = "/home/tlangloi/work/sandbox/rigid-ipc-chrono/"
    # out_folder = "xyz"
    # json_file = "xyz.json"
    out_folder = "slopeTest_highSchoolPhysics_mu=0.5"
    json_file = "fixtures/3D/friction/incline-plane/slopeTest_highSchoolPhysics_mu=0.5.json"
    # out_folder = "large-mass-ratio"
    # json_file = "fixtures/3D/unit-tests/large-mass-ratio.json"
    # out_folder = "spike-and-wedge"
    # json_file = "fixtures/3D/unit-tests/erleben/spike-and-wedge.json"
    # out_folder = "wedges"
    # json_file = "fixtures/3D/unit-tests/erleben/wedges.json"
    # out_folder = "vertex-face"
    # json_file = "fixtures/3D/unit-tests/vertex-face.json"
    # out_folder = "arch-101-stones"
    # json_file = "fixtures/3D/friction/arch/arch-101-stones.json"
    # out_folder = "edge-edge"
    # json_file = "fixtures/3D/unit-tests/edge-edge.json"
    # out_folder = "spikes"
    # json_file = "fixtures/3D/unit-tests/erleben/spikes.json"
    # out_folder = "internal-edges"
    # json_file = "fixtures/3D/unit-tests/erleben/internal-edges.json"
    # out_folder = "wedge-in-crack"
    # json_file = "fixtures/3D/unit-tests/erleben/wedge-in-crack.json"
    # out_folder = "spike-in-crack"
    # json_file = "fixtures/3D/unit-tests/erleben/spike-in-crack.json"
    # out_folder = "5-cubes"
    # json_file = "fixtures/3D/unit-tests/5-cubes.json"
    # out_folder = "screw"
    # json_file = "fixtures/3D/mechanisms/screw.json"
    # out_folder = "chain-net-8x8"
    # json_file = "fixtures/3D/chain/chain-net/chain-net-8x8.json"
    # out_folder = "chain-net-4x4"
    # json_file = "fixtures/3D/chain/chain-net/chain-net-4x4.json"
    # out_folder = "10-links"
    # json_file = "fixtures/3D/chain/10-links.json"

    # The path to the Chrono data directory containing various assets (meshes, textures, data files)
    # is automatically set, relative to the default location of this demo.
    # If running from a different directory, you must change the path to the data directory with:
    chrono.SetChronoDataPath(root_path)
    mesh_path = os.path.join(root_path, "fixing-collisions/meshes")
    #mesh_path = ""

    out_path = os.path.join(root_path, out_folder)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    with open(os.path.join(root_path, "fixing-collisions", json_file)) as f:
        data = json.load(f)

    system = chrono.ChSystemNSC()
    if "gravity" in data["rigid_body_problem"]:
        g = data["rigid_body_problem"]["gravity"]
        system.Set_G_acc(chrono.ChVectorD(g[0], g[1], g[2]))

    # system.SetSolverType(chrono.ChSolver.Type_APGD)
    system.SetSolverType(chrono.ChSolver.Type_BARZILAIBORWEIN)
    system.SetSolverMaxIterations(15000)

    # Set the global collision margins. This is expecially important for very large or
    # very small objects. Set this before creating shapes. Not before creating system.
    # chrono.ChCollisionModel.SetDefaultSuggestedEnvelope(0.001)
    # chrono.ChCollisionModel.SetDefaultSuggestedMargin(0.001)
    chrono.ChCollisionModel.SetDefaultSuggestedEnvelope(0.1)
    chrono.ChCollisionModel.SetDefaultSuggestedMargin(0.00001)

    # ---------------------------------------------------------------------
    #
    #  Create the simulation system and add items
    #

    # Create a contact material (with default properties, shared by all collision shapes)
    contact_material = chrono.ChMaterialSurfaceNSC()

    has_friction = False
    if "coefficient_friction" in data["rigid_body_problem"]:
        contact_material.SetFriction(
            data["rigid_body_problem"]["coefficient_friction"])
        has_friction = True

    meshes = []

    for body in data["rigid_body_problem"]["rigid_bodies"]:
        mpath = os.path.join(mesh_path, body["mesh"])

        if body["mesh"] == "plane.obj":
            pos = body["position"]
            dim = body["dimensions"] if "dimensions" in body else [
                10, 0.005, 10]
            rot = body["rotation"] if "rotation" in body else [0.0, 0.0, 0.0]
            rot = Rotation.from_euler('xyz', rot, degrees=True)
            rot = rot.as_quat()

            floor = chrono.ChBodyEasyBox(
                dim[0], dim[1], dim[2], 1, True, True, contact_material)
            floor.SetPos(chrono.ChVectorD(pos[0], pos[1], pos[2]))
            floor.SetRot(chrono.ChQuaternionD(rot[3], rot[0], rot[1], rot[2]))
            if "is_dof_fixed" in body and body["is_dof_fixed"]:
                floor.SetBodyFixed(True)
            system.Add(floor)
            meshes.append((np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int)))
            continue

        if "scale" in body and body["scale"] != 1:
            v, f = igl.read_triangle_mesh(mpath)
            v *= body["scale"]
            igl.write_triangle_mesh("tmp.obj", v, f)
            mpath = "tmp.obj"

        density = 1000
        if "density" in body:
            density = body["density"]

        cbody = chrono.ChBodyEasyMesh(mpath,  # mesh filename
                                      density,             # density kg/m^3
                                      True,             # automatically compute mass and inertia
                                      True,             # visualize?>
                                      True,             # collide?
                                      contact_material,  # contact material
                                      )

        cbody.Update()
        V, F = igl.read_triangle_mesh(mpath)
        t, R = get_trafos(cbody)
        print(t)
        print(R)
        V = V - t
        V = V @ R
        meshes.append((V, F))

        #rot = Rotation.from_euler('zyx', body["rotation"], degrees=True)
        rot2 = Rotation.from_matrix(R)
        rotTmp = body["rotation"] if "rotation" in body else [0.0, 0.0, 0.0]
        rot = Rotation.from_euler('xyz', rotTmp, degrees=True) * rot2
        rot = rot.as_quat()
        cbody.SetRot(chrono.ChQuaternionD(rot[3], rot[0], rot[1], rot[2]))

        if "position" in body:
            pos = body["position"]
            pos = [x - y for (x, y) in zip(pos, t)]
            cbody.SetPos(chrono.ChVectorD(pos[0], pos[1], pos[2]))

        # print(body)

        # Need to fix coordinate systems
        # if "linear_velocity" in body:
        #     vel = body["linear_velocity"]
        #     cbody.SetPos_dt(chrono.ChVectorD(vel[0], vel[1], vel[2]))

        # if "angular_velocity" in body:
        #     vel = body["angular_velocity"]
        #     cbody.SetWVel_loc(chrono.ChVectorD(vel[0], vel[1], vel[2]))

        is_fixed = False
        if "type" in body:
            is_fixed = body["type"] == "static"
        elif "is_dof_fixed" in body:
            if isinstance(body["is_dof_fixed"], list):
                for k in range(min(6, len(body["is_dof_fixed"]))):
                    if not body["is_dof_fixed"][k]:
                        print("warning only all dofs fixed supported")
                    else:
                        is_fixed = True
            else:
                is_fixed = body["is_dof_fixed"]
        cbody.SetBodyFixed(is_fixed)

        system.Add(cbody)

    application = chronoirr.ChIrrApp(
        system, 'runner', chronoirr.dimension2du(1024, 768))

    application.AddTypicalCamera(chronoirr.vector3df(
        0.5, 0.5, 1), chronoirr.vector3df(0, 0, 0))
    application.AddTypicalLights()
    application.AssetBindAll()
    application.AssetUpdateAll()

    application.SetTimestep(data["timestep"])
    tt = 0
    index = 0

    save_mesh(application, out_path, 0)

    while(application.GetDevice().run()):
        application.BeginScene()
        application.DrawAll()
        application.DoStep()
        # if index == 0:
        #     application.DoStep()
        #     index += 1
        application.EndScene()

        tt += application.GetTimestep()
        index += 1
        print("t = {}".format(tt))
        save_mesh(application, out_path, index)

        if tt >= data["max_time"]:
            break
