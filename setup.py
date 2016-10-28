#!/usr/bin/env python

from setuptools import setup, find_packages
#from distutils.core import setup

setup(name='SimulationPipeline',
      version='1.1',
      description='',
      author='Kasper Munch',
      author_email='kaspermunch@birc.au.dk',
      #packages=['SimulationPipeline'],
      packages=find_packages(),
	  entry_points = { 'console_scripts': [ 'run_simulation_pipeline = SimulationPipeline.Utils:runSimulationsWithCoaSimAndILS09Script',
                                            'run_ils16_simulation_pipeline = SimulationPipeline.Utils:runSimulationsWithCoaSimAndILS16Script',
                                            'run_simulation_pipeline_with_ctmc = SimulationPipeline.Utils:runSimulationsWithCoaSimAndILSCTMCScript', ]
           }
     )
