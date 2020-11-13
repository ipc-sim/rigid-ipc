import pathlib

fixture_dir = pathlib.Path(__file__).resolve().parents[2] / "fixtures"

scenes = [
    "3D/5-cubes.json",
    "3D/simple",
    "3D/erleben",
    # "3D/armadillos.json",
    "3D/bat.json",
    # "3D/blender/blender-1000.json",
    "3D/chain",
    "3D/compactor",
    "3D/rotation",
    "3D/dzhanibekov.json",
    "3D/edges",
    "3D/example.json",
    "3D/fracture",
    "3D/gears",
    "3D/large-mass-ratio.json",
    "3D/lever-arm.json",
    "3D/pendulum",
    # "3D/piles/cone-armadillos.json",
    "3D/piles/cone.json",
    "3D/piles/plane.json",
    "3D/tessellated-plane",
    "3D/tet-corner.json",
    "3D/tunnel",
    "3D/tunneling.json",
    "3D/turntable/turntable-mu=0_0.json",
    "3D/piles/cone-bunnies-lowpoly.json",
    "3D/piles/cone-bunnies.json",
    "3D/blender/blender-1.json",
    "3D/blender/blender-10.json",
    "3D/blender/blender-25.json",
    "3D/blender/blender-100.json",
    "3D/blender/blender-200.json",
    "3D/screw.json",
    "3D/wrecking-ball.json",
]

print("-i", end=" ")

for scene in scenes:
    if(scene != ""):
        print(fixture_dir / (scene), end=" ")

print("-o no-friction.csv", end="")
