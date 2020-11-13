import pathlib

fixture_dir = pathlib.Path(__file__).resolve().parents[2] / "fixtures"

scenes = [
    "3D/friction",
    "3D/rolling",
    "3D/turntable/turntable-mu=0_1.json",
    "3D/turntable/turntable-mu=0_5.json",
]

print("-i", end=" ")

for scene in scenes:
    if(scene != ""):
        print(fixture_dir / scene, end=" ")

print("-o friction.csv", end="")
