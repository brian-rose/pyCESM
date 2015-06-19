from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

import pyCESM

setup(name='pyCESM',
      version=pyCESM.__version__,
      description='Package for interactive work with the Community Earth System Model',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
      ],
      keywords='climate modeling modelling model gcm',
      url='http://github.com/brian-rose/pyCESM',
      author='Brian E. J. Rose',
      author_email='brose@albany.edu',
      license='MIT',
      packages=['pyCESM'],
      install_requires=[
          'numpy',
      ],
      include_package_data=True,
      zip_safe=False)
