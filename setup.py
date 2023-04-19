from setuptools import setup, find_packages

# Load the version number from inside the package
exec(open('workflows/__version__.py').read())

# Use the README as the long_description
with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='feature-workflows',
    version='',
    packages=find_packages(),
    url='https://github.com/BSDExabio/feature-workflows',
    license='',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Russell B. Davidson, Mark Coletti, Ada Sedova',
    author_email='davidsonrb@ornl.gov, colettima@ornl.gov, sedovaaa@ornl.gov',
    description='Protein feature extraction workflows',
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'dask[complete]',  # Used for parallel and distributed algorithms
        'distributed',  # Used for parallel and distributed algorithms
        'matplotlib',  # Used in visualizations
        'numpy',  # Used for vector math
        'pandas',  # Used to process CSV output
        'rich',  # Used for pretty printing logs etc.
        'scipy',
        'toolz'  # All round useful extension to itertools and functools
    ]
)
