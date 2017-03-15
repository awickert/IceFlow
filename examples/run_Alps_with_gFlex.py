import iceflow
ic = iceflow.IceFlow()
self = ic # For interactive troubleshooting

import numpy as np
import scipy
from scipy.io import loadmat
from scipy.sparse import linalg
import os
import sys
from matplotlib import pyplot as plt

ic.useGRASS = True
ic.location='AlpenJuergen'
# CHANGE GISBASE FOR YOUR COMPUTER!
ic.gisbase = '/usr/local/src/grass7_trunk/dist.x86_64-unknown-linux-gnu'

ic.run_length_years = 200. # years
ic.t_start_years = 0. # years
ic.dt_years = 1
ic.dt_is_variable = True # if this is set, the above dt_years still counts
                         # for the first time step (but there is a default
                         # that will take over if I unset it)
ic.dt_max = 5 # maximum permissible timestep

ic.elevation = 'BedrockJuergenAlgorithm'

ic.north = 5400000
ic.south = 4820000
ic.west = 170000
ic.east = 1100000
ic.dx = 2500
ic.dy = 2500

# Set this up to automatically number IceFlow outputs using glob
ic.output_filename=None
ic.output_figure=None
ic.plot_at_end_flag=False
ic.plot_during_run_flag = True
#ic.plot_t_years = ic.run_length_years
ic.boundary_condition = 'Dirichlet0'

ic.GRASS_raster_ice_extent = 'IceExtentAlpenLGM'
#ic.GRASS_raster_ice_extent = 'LGM_Alpen_noholes'

ic.verbose = True

ic.record_frequency_years = 100.

# Flexure
ic.isostatic = True
#ic.ElasticThickness = 20000
ic.ElasticThickness = 'EETEurope_meters' # Cell size is 20000 m; this is
                                         # minimum grid size for flexure
# Bounding box for flexure
ic.flex.Quiet=True
ic.flex.north = 5600000
ic.flex.south = 4000000
ic.flex.west = -300000
ic.flex.east = 1500000
ic.flex.dx = 50000
ic.flex.dy = 50000

ic.isostatic_response_time_scale = 100. # Intentionally short to come to 
                                        # equilibrium quickly

# In order to loop over possible mass balances, just create a "for" loop with 
# these, and then output the resultant values

ic.mass_balance_parameterization = 'ELA'
ic.ELA = 2300
ic.dbdz_per_year = 1E-3 # 1 m/yr over 1 km
ic.b_maximum_per_year = 0.3 # [m/yr]

ic.initialize()
ic.run()
ic.finalize()

# How well the ice model fits the data is determined by these

