import numpy
import pathlib

import rigidipc


rigidipc.set_logger_level(rigidipc.LoggerLevel.error)

sim = rigidipc.Simulation()

root_path = pathlib.Path(__file__).parents[1]
scene = root_path / "fixtures" / "3D" / "unit-tests" / "5-cubes.json"
sim.load_scene(str(scene))
print(f"timestep={sim.timestep:g}")
sim.max_simulation_steps = 50

# Make body 0 a kinematic object
body = sim.bodies()[0]
body.name = "kinematic_cube"
body.type = rigidipc.KINEMATIC
# Prescribe the kinematic motion of the body for the first 50 steps
body.kinematic_poses = [
    rigidipc.Pose(numpy.array([0, i / 12.5, 0]), numpy.array([i / 12.5, 0, 0]))
    for i in range(50)]

# Run the simulation
output_path = pathlib.Path(__file__).parent / "output"
rigidipc.set_profiler_output_directory(str(output_path))
sim.run(str(output_path / "sim.json"))

for body in sim.bodies():
    print(f"{body.name}: pose={body.pose}, group_id={body.group_id}, type={body.type}")
