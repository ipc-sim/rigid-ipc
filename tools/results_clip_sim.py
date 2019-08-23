"""Remove the simulation frames after n."""

import sys
import json

filename = sys.argv[1]

with open(filename) as file:
    sim = json.load(file)

n = int(sys.argv[2])
sim["animation"]["vertices_sequence"] = (
    sim["animation"]["vertices_sequence"][:n])
sim["animation"]["state_sequence"] = sim["animation"]["state_sequence"][:n]

with open("clipped_sim.json", "w") as file:
    json.dump(sim, file)
