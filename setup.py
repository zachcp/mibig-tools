from setuptools import setup, find_packages

setup(
    name='mibig-tools',
    version='0.1.13',
    install_requires=[
        'Click >= 4.0',
        'Biopython >=1.6.5',
        'toolz >= 0.7.2',
        'pandas >=0.16'
    ],
    packages=find_packages(),
    entry_points='''
        [console_scripts]
        mibig_process_folder = mibig_tools.process_all_mibig:process_mibig_cluster_folder
    ''',
)
