from setuptools import setup, find_packages

setup(
    name='cleanpipe',
    version='0.1',
    packages=find_packages(),
    description='Functions that make computational biology pipelines easyer to build and understand.',
    author='Henrique Rigitano',
    author_email='henrique.rigitano@alumni.usp.br',
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
        'py3Dmol',
        'scipy',
        'scikit-learn',
    ],

    include_package_data=True,  # this triggers MANIFEST.in to handle recursive inclusion of non-python files


)
