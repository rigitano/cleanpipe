from setuptools import setup, find_packages

setup(
    name='cleanpipe',
    version='0.1',
    packages=find_packages(),
    description='A custom package for clean code',
    author='Henrique Rigitano',
    author_email='henrique.rigitano@gmail.com',
    url='https://github.com/rigitano/cleanpipe',
    install_requires=[
        'PeptideBuilder', 
        'Bio',
        'Geometry',
        'mdtraj',
        'pandas',
        'numpy',
        'requests',
        'matplotlib',
        'seaborn',
        'MDAnalysis',
    ],
    include_package_data=True, # Normally pip install consider only .py files. this will alow it to see the mdp files in the mdp folder
    package_data={
        'cleanpipe.mdp': ['*.mdp']
    },
)