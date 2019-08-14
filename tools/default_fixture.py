"""Function to get the default dictionary for fixtures."""


def generate_default_fixture():
    """Create the default fixture as a dictionary."""
    return {
        "scene_type": "distance_barrier_rb_problem",
        "collision_solver": "barrier_solver",
        "timestep_size": 1e-2,
        "distance_barrier_constraint": {
            "custom_initial_epsilon": 0.1,
            "detection_method": "hash_grid",
            "use_distance_hashgrid": True,
            "active_constraint_scale": 1.01,
            "custom_hashgrid_cellsize": -1
        },
        "barrier_solver": {
            "min_barrier_epsilon": 1e-2
        },
        "rigid_body_problem": {
            "gravity": [0.0, 0.0, 0.0],
            "coefficient_restitution": 0.0,
            "rigid_bodies": [],
        }}
