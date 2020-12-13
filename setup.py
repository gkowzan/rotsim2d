from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(
    name='rotsim2d',
    version='0.0.1',
    description="Simulation of rotationally-resolved rovibrational 2D spectra",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gitlab.com/allisonlab/mdcs/rotsim2d',
    author='Grzegorz Kowzan',
    author_email='grzegorz.kowzan@stonybrook.edu',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Physics'],
    python_requires='>=3.5',
    install_requires=['numpy', 'scipy', 'anytree',
                      'pywigxjpf @ git+ssh://git@gitlab.com:allisonlab/mdcs/pywigxjpf.git@master',
                      'shed @ git+ssh://git@gitlab.com:allisonlab/mdcs/shed.git@master',
                      'spectroscopy @ git+ssh://git@gitlab.com:allisonlab/mdcs/spectroscopy.git@master'],
    extras_require={'doc': ['sphinx', 'sphinx_rtd_theme', 'numpydoc']},
    packages=find_packages(),
)
