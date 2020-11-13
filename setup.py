from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='emptydrops',
    version='0.0.5',
    author='nh3',
    author_email='nh3@users.noreply.github.com',
    description='Python implementation of emptydrops() in CellRanger v3.0.2',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/nh3/emptydrops',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'h5py',
        'lz4',
        'pandas',
        'scipy',
    ],
)
