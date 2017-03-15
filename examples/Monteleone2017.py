import iceflow
ic = iceflow.IceFlow()
self = ic

import numpy as np
import scipy
from scipy.io import loadmat
from scipy.sparse import linalg
import os
import sys
from matplotlib import pyplot as plt

ic.useGRASS = True
ic.location='Monteleone2017_projected'
ic.gisbase = '/usr/local/src/grass7_trunk/dist.x86_64-unknown-linux-gnu'

ic.run_length_years = 5000. # years
ic.t_start_years = 0. # years
ic.dt_years = 0.05
ic.record_frequency_years = 50

ic.elevation = 'gebco_08'

ic.north = 6700000
ic.south = 6040000
ic.west = 400000
ic.east = 790000
ic.dx = 1000
ic.dy = 1000

ic.mass_balance_parameterization = 'ELA'
ic.ELA = 1200
ic.dbdz_per_year = 1E-3 # 1 m/yr over 1 km -- EDIT THIS TO FIT DATA!
ic.b_maximum_per_year = .3 # [m/yr]

# Set this up to automatically number IceFlow outputs using glob
ic.output_filename=None
ic.output_figure=None
ic.plot_at_end_flag=False
ic.plot_during_run_flag = True
ic.OutNameGRASS = 'ELA'+str(ic.ELA)+'_tmax'+str(ic.run_length_years)+'_dt'+str(ic.dt_years)+'_dx'+str(ic.dx)
#ic.plot_t_years = ic.run_length_years
ic.boundary_condition = 'Dirichlet0'

ic.GRASS_raster_ice_extent = 'FakeMeasuredExtents'

ic.verbose = False

# Flexure
ic.isostatic = False

ic.initialize()
ic.run()
ic.finalize()


