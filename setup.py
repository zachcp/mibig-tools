from setuptools import setup, find_packages

setup(
    name='mibig-tools',
    version='0.1.0',
    install_requires=[
        'Click >= 0.4.0',
        'Biopython >=1.6.5'
    ],
    packages=find_packages(),
    entry_points='''
        [console_scripts]
        mibig_process_cluster = mibig_tools.process_cluster:general_cluster_data
    ''',
)
