
import ez_setup
ez_setup.use_setuptools()

import sys, glob
from setuptools import setup, Extension, find_packages

if sys.version < '2.7':
    print "This package requires Python version 2.7 or higher. Aborting."
    print "Visit http://www.python.org/download for a more recent version of Python."
    print "Aborting."
    sys.exit()

version = "1.0.0"
package_name = "SimulationPipeline"

setup(name=package_name,
      test_suite='tests',
      version=version,
      description='',
      long_description='',
      author='Kasper Munch',
      author_email='kaspermunch@birc.au.dk',
      url='http://users-cs.au.dk/kmt',
      packages = find_packages(exclude=['ez_setup']),
      package_dir = {package_name: package_name},
      include_package_data = True,
      entry_points = { 'console_scripts': [ 'run_simulation_pipeline = SimulationPipeline.Utils:runSimulationsWithCoaSimScript', ],
                       },
      )
