#!/bin/bash
python ./tools/generate_chainmail_fixture.py 10 \
    ./comparisons/Box2D/example.box2d.json

./build/release/FixingCollisions_ngui ./comparisons/Box2D/example.box2d.json \
    ./results/comparisons/Box2D/chainmail/ours --num-iterations 1000
python ./tools/results_to_vtk_files.py \
    ./results/comparisons/Box2D/chainmail/ours/sim.json \
    ./results/comparisons/Box2D/chainmail/ours

./comparisons/Box2D/run.sh
mv ./results/comparisons/Box2D/results.json \
    ./results/comparisons/Box2D/chainmail/box2d/sim.json
python ./tools/results_to_vtk_files.py \
    ./results/comparisons/Box2D/chainmail/box2d/sim.json \
    ./results/comparisons/Box2D/chainmail/box2d
