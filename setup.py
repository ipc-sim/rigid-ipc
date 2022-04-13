import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements\n"
        " in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

requirements = ["autolab-core", "numpy"]
exec(open("src/rigidipc/version.py").read())

from setuptools import find_packages

setup(
    name="rigidipc",
    version=__version__,
    # description="PhoXi python driver",
    # long_description="Python bindings for Photoneo PhoXi camera driver",
    # author="Mike Danielczuk",
    # author_email="michael.danielczuk@gmail.com",
    # license="MIT",
    install_requires=requirements,
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_install_dir="src/rigidipc",
    include_package_data=True,
    python_requires=">=3.6",
)
