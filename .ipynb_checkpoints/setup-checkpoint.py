# setup.py

from setuptools import setup, find_packages

with open("README.md", "r") as f:
    description = f.read()


setup(
    name='cfdpre',
    version='0.1',
    description='CFD PreProcessing Library',
    author='Pushkar Sheth',
    author_email='siglyserdev@gmail.com',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'CoolProp'

        
    ],

    long_description=description,
    long_description_content_type="text/markdown",
)