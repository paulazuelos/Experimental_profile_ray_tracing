#!usr/bin/python
# -*- coding: utf-8 -*-

"""Package installation setup !"""

import os
import subprocess
from setuptools import setup, find_packages

version = '0.0.1'
sha = 'Unknown'
package_name = 'ray_tracing_lib'

cwd = os.path.dirname(os.path.abspath(__file__))

print("Building wheel {}-{}".format(package_name, version))


def write_version_file():
    version_path = os.path.join(cwd, 'ray_tracing_lib', 'version.py')
    with open(version_path, 'w') as f:
        f.write("__version__ = '{}'\n".format(version))


write_version_file()

with open('README.md') as f:
    readme = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name=package_name,
    version=version,
    author='Paul Azuelos',
    description='Ray tracing from Experiemental profile',
    long_description=readme,
    long_description_content_type="text/markdown",
    url='',
    download_url='',
    license='Copyright',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Copyright License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    keywords=['optics', 'photonics', 'lens', 'simulation'],
    packages=find_packages(exclude=['test', 'docs']),
    zip_safe=True,
    python_requires='>=3.6.0',
    include_package_data=True,
    install_requires=requirements,
    package_data={'': ['LICENSE']},
)