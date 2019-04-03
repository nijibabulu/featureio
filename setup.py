import io
from setuptools import setup


setup(
    name='featureio',
    version='0.0.1',
    url='',
    license='MIT',
    author='Bob Zimmermann',
    author_email='robert.zimmermann@univie.ac.at',
    description='A genomic feature format parser and converter',
    long_description=io.open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    packages=["featureio"],
    test_require=["pytest"]
)


