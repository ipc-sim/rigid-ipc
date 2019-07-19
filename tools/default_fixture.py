"""Function to get the default dictionary for fixtures."""


def generate_default_fixture():
    """Create the default fixture as a dictionary."""
    return {
        "scene_type": "rigid_body_problem",
        "collision_solver": "barrier_solver",
        "rigid_body_problem": {
            "rigid_bodies": [],
            "constraint": "distance_barrier_constraint",
            "update_constraint_set": True,
            "use_chain_functional": False,
            "gravity": [0.0, 0.0, 0.0],
            "collision_eps": 2.0
        },
        "particles_problem": {
            "vertices": [],
            "edges": [],
            "velocities": [],
            "x_fixed": [],
            "y_fixed": [],
            "use_mass_matrix": True,
            "constraint": "distance_barrier_constraint",
            "update_constraint_set": True,
            "gravity": [0.0, 0.0],
            "collision_eps": 2.0
        },
        "barrier_solver": {
            "inner_solver": "newton_solver",
            "min_barrier_epsilon": 1e-5,
            "max_iterations": 0
        },
        "gradient_descent_solver": {
            "absolute_tolerance": 1e-5,
            "min_step_length": 1e-12,
            "max_iterations": 3000
        },
        "newton_solver": {
            "absolute_tolerance": 1e-5,
            "min_step_length": 1e-12,
            "max_iterations": 3000
        },
        "bfgs_solver": {
            "absolute_tolerance": 1e-5,
            "min_step_length": 1e-12,
            "max_iterations": 3000
        },
        "time_barrier_constraint": {
            "initial_epsilon": "min_toi",
            "custom_initial_epsilon": 1.0,
            "detection_method": "hash_grid",
            "extend_collision_set": True
        },
        "distance_barrier_constraint": {
            "custom_initial_epsilon": 0.5,
            "detection_method": "hash_grid",
            "extend_collision_set": False
        },
        "timestep_size": 0.1,
        "viewport_bbox": {"min": [0, 0], "max": [0, 0]}
    }
