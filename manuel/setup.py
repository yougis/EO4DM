from setuptools import setup, find_packages
import os
import re

directory = os.path.dirname(os.path.abspath(__file__))


# Extract version information
path = os.path.join(directory, 'dmpipeline', '__init__.py')
with open(path) as read_file:
    text = read_file.read()
pattern = re.compile(r"^__version__ = ['\"]([^'\"]*)['\"]", re.MULTILINE)
version = pattern.search(text).group(1)

# Extract long_description
path = os.path.join(directory, 'README.md')
with open(path) as read_file:
    long_description = read_file.read()


setup(
    name='dmpipeline',
    version=version,
    author="Mathis Neuhauser, Thomas Tilak, OEIL",
    description="Pipeline to compute Drought Monitoring Indicators",
    long_description=long_description,
    packages=find_packages(include=['dmpipeline',
                                    'dmpipeline.GEE_Processing',
                                    'dmpipeline.GEE_Processing.gee_accounts',
                                    'dmpipeline.DROUGHT_Processing',
                                    'dmpipeline.ALERT_Processing',
                                    'dmpipeline.ALERT_Processing.ftp_accounts',
                                    'dmpipeline.GEOSTATS_Processing']),
    package_data={
        'dmpipeline.GEE_Processing.gee_accounts':['*.json'],
        'dmpipeline.ALERT_Processing.ftp_accounts':['*.json']
        },
    # install_requires=[
    #     'rasterio',
    #     'fiona',
    #     'earthengine-api'
    # ],
    entry_points={
        'console_scripts': [
            'globalDM = dmpipeline.GlobalDrought_ProcessingChain:main',
            'localDM = dmpipeline.LocalDrought_ProcessingChain:main',
            'alertDM = dmpipeline.AlertDrought_ProcessingChain:main',
        ]
    }
)

