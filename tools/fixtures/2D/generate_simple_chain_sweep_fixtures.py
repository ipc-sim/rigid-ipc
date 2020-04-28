#!/usr/local/bin/python3
"""
Script to generate a fixture of a chain with N simple links.

Usage: python generate_simple_chainmail_fixture.py N
"""

import json
import pathlib
import sys

import numpy as np

from fixture_utils import *


def link(a, b, c, d, e):
    #           e
    #         _____
    #       b|     |c
    # a|_____|  |__'__
    #  |  d  |  |  ,
    #        |_____|
    #            f

    V = np.array([[0.0, -a], [0.0, 0.0], [0.0, a], [d, -b], [d, 0.0], [d, b],
                  [d + e, -b], [d + e, -b + c], [d + e, +b - c], [d + e, +b]],
                 dtype=np.float64)
    E = np.array([[0, 1], [1, 2], [3, 4], [4, 5], [6, 7], [8, 9], [1, 4],
                  [3, 6], [5, 9]],
                 dtype=np.int32)

    return V, E

def generate_fixture(n_links: int, d: float, s: float) -> dict:
    a = 0.4
    b = 0.5
    c = 0.4
    # d = 0.5
    e = 0.5
    f = 0.5

    vertices, edges = link(a, b, c, d, e)
    dx = np.array([d + e - f * d, 0.0], dtype=np.float64)
    vertices[:,1] *= s



    """Generate a fixture of a chain with N simple links."""
    fixture = generate_default_fixture()
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]
    fixture["rigid_body_problem"]["gravity"] = [0.0, -9.8, 0.0]
    fixture["rigid_body_problem"]["coefficient_restitution"] = 0
    for i in range(n_links):
        rigid_bodies.append({
            "vertices": (vertices + dx * i).tolist(),
            "edges": edges.tolist(),
            "velocity": [0.0, 0.0 if i else 0.0, 0.0],
            "is_dof_fixed": [i == 0, i == 0, i == 0],
        })
    return fixture


def main() -> None:
    """Parse command-line arguments to generate the desired fixture."""
    assert len(sys.argv) >= 2

    if len(sys.argv) > 2:
        out_path = pathlib.Path(sys.argv[2])
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        directory = pathlib.Path(
            __file__).resolve().parents[1] / "fixtures" / "chain_sweep"
        directory.mkdir(parents=True, exist_ok=True)
        

    
    n_links = int(sys.argv[1])
    for i, d in enumerate([0.5, 0.2, 0.1]):
        for j, s in enumerate([1., 0.7, 0.4]):
            fixture = generate_fixture(n_links, d, s)
            out_path = directory / f"simple_{n_links:d}_link_chain_x_{i}_y_{j}.json"
            save_fixture(fixture, out_path)


if __name__ == "__main__":
    main()

