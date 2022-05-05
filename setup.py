"""
Apply MESMER segmentation to NanoString SMI CellComposition images
"""
from setuptools import find_packages, setup
from distutils.core import setup, Extension

dependencies = ['click', 'numpy<1.20.0', 'scikit-image<0.19.0', 'matplotlib',  'deepcell', 'PyYAML']
 
setup(
    name='mesmer_segment',
    version='0.1.0',
    url='https://github.com/Jacobog02/mesmer_segment]',
    license='MIT',
    author='Jacob Gutierrez',
    author_email='jacobog@stanford.edu',
    description='apply mesmer segmentation to SMI CellComposite images',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    package_data={},        ## Point to extra package data
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'mesmer_segment = mesmer_segment.cli:main',
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
         'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
        # 'Development Status :: 3 - Alpha',
        # 'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
    
)

