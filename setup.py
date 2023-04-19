from setuptools import setup

# Load the version number from inside the package
exec(open('workflows/__version__.py').read())

# Use the README as the long_description
with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='feature-workflows',
    version='',
    packages=[''],
    url='',
    license='',
    long_description=long_description,
    author='Russell B. Davidson, Mark Coletti, Ada Sedova',
    author_email='davidsonrb@ornl.gov, colettima@ornl.gov, sedovaaa@ornl.gov',
    description='Protein feature extraction workflows'
)
