import pathlib


def print_fixtures(name, scenes):
    fixture_dir = (
        pathlib.Path(__file__).resolve().parents[2] / "fixtures" / "3D")

    print("-i", end=" ")

    for scene in scenes:
        for p in fixture_dir.glob(scene + ".json"):
            if p.exists():
                print(p, end=" ")
            else:
                raise Exception(f"{p} does not exist!")

    print(f"-o {name}.csv", end="")
