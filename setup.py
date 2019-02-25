#!/usr/bin/env python2

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FasTosh",
    version="0.1",
    author="Pauline Barbet",
    author_email="arnaud.felten@anses.fr",
    description="FasTosh: Fast Taxonomy with mash",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/afelten-Anses/QuickPhylo",
    packages=setuptools.find_packages(),
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Programming Language :: Python :: 2",
        "Operating System :: POSIX :: Linux",
    ],
        scripts=["FasTosh",
             "comparaison.py",
             "tree.py"
             ],
    include_package_data=True,
    install_requires=['dendropy',
                      'mash',
                      ], 
    zip_safe=False,

)
