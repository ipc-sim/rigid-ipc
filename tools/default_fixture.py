"""Function to get the default dictionary for fixtures."""


def generate_default_fixture():
    """Create the default fixture as a dictionary."""
    return {
        "scene_type": "rigid_body_problem",
        "rigid_body_problem": {
            "rigid_bodies": [],
            "gravity": [0.0, 0.0, 0.0],
        }}
