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
        'nglview',
    ],

    include_package_data=True,  # this triggers MANIFEST.in handling

    # not needed because MANIFEST.in handles recursive inclusion. this is just a reminder that doing this is necessary, otherwise just py files will be considered
    package_data={
        "cleanpipe": [
            "mdp/*",
            "bash/*",
            "tcl/*",
            "USEFUL_SOLVENTS/*",
            "USEFUL_MOLECULES/*",
            "USEFUL_FORCEFIELDS/*",
        ]
    },
)
