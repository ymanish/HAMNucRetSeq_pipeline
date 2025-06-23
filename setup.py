from setuptools import setup, find_packages

setup(
        name="hammoud_reanalysis",
        version='0.1',
        packages=find_packages(),
        author='Manish Yadav',
        author_email='manish20072013 at gmail dot com',
        description='A pipeline for processing the raw sequencing data of retented nucleosomes in Sperm cells from Hmmoud et al. 2009',
        long_description=open('README.md', encoding='utf-8').read(),
        long_description_content_type='text/markdown',
        license='GPL3', 
        python_requires=">=3.10",
        
    )