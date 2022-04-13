import numpy
import pathlib

import rigidipc as ri

ri.set_logger_level(ri.LoggerLevel.error)

sim = ri.Simulation()

scene = "fixtures/3D/unit-tests/5-cubes.json"
sim.load_scene(str(scene))
print(f"timestep={sim.timestep:g}")
sim.max_simulation_steps = 50

# Make body 0 a kinematic object
body = sim.bodies()[0]
body.name = "kinematic_cube"
body.type = ri.KINEMATIC
# Prescribe the kinematic motion of the body for the first 50 steps
body.kinematic_poses = [
    ri.Pose(numpy.array([0, i / 12.5, 0]), numpy.array([i / 12.5, 0, 0]))
    for i in range(50)]

# Run the simulation
output_path = "output"
ri.set_profiler_output_directory(str(output_path))
sim.run(f"{output_path}/sim.json")

for body in sim.bodies():
    print(f"{body.name}: pose={body.pose}, group_id={body.group_id}, type={body.type}")
