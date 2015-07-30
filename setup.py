#!/usr/bin/env python

from distutils.core import setup

setup(name='SimulationPipeline',
      version='1.0',
      description='',
      author='Kasper Munch',
      author_email='kaspermunch@birc.au.dk',
      packages=['SimulationPipeline'],
	  entry_points = { 'console_scripts': [ 'run_simulation_pipeline = SimulationPipeline.Utils:runSimulationsWithCoaSimScript',
                                            'run_simulation_pipeline_with_ctmc = SimulationPipeline.Utils:runSimulationsWithCoaSimAndILSCTMCScript', ]
           }
     )
