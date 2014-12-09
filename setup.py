import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'endseq', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
endseq - analysis pipeline for poly(A)-tail length and 3`end analysis
"""

setup(
        name="endseq",
        version=version,
        install_requires=['matplotlib >1.3.0',
                          'pandas >0.14.0',
                          'numpy > 1.8.0',
                          'HTSeq',
                          'seaborn',
                          'brewer2mpl'],
        requires = ['python (>=2.7, <3.0)'],
        packages=['endseq',
                  'endseq.scripts'],
        author="Mohan Bolisetty",
        description='A toolset for working with endseq data',
        long_description=long_description,
        url="",
        package_dir = {'endseq': "endseq"},
        package_data = {'endseq': []},
        zip_safe = False,
        include_package_data=True,
##        scripts = ['endseq/scripts/endseq-script'],
        entry_points = {
            'console_scripts' : [
                 'endseq = endseq.endseq_main:main', 
            ],
        },  
        author_email="mohanbolisetty@gmail.com",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
