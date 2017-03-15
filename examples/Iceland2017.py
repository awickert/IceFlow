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
ic.location='Iceland'
ic.gisbase = '/usr/local/src/grass7_trunk/dist.x86_64-unknown-linux-gnu'

ic.run_length_years = 1000. # years
ic.t_start_years = 0. # years
ic.dt_years = 1
ic.record_frequency_years = 50

ic.elevation = 'elev'

ic.north = 792000
ic.south = 240000
ic.west = 108000
ic.east = 868000
ic.dx = 4000
ic.dy = 4000

"""
ic.mass_balance_parameterization = 'TP_PDD'
# Goes to Ta variable
ic.temperature_melt_season = 'TwentiethCenturyTemp_NCAR_CIRES_JJA_degC'
ic.melt_season_length_days = 120
#self.Pa /= 1000. / self.secyr # e.g., PRISM is in mm/yr, change to m/s.
# CONVERT THIS IN GRASS BEFOREHAND!!! <----------------------------------------
# Goes to Pair variable
ic.precipitation_solid = 'TwentiethCenturyPrecip_NCAR_CIRES_annual'
ic.T_correction = -5
ic.P_factor = 1
ic.b_maximum_per_year = 1 # [m/yr]
"""

ic.mass_balance_parameterization = 'ELA'
ic.ELA = 800
ic.dbdz_per_year = 1E-3 # 1 m/yr over 1 km -- EDIT THIS TO FIT DATA!
ic.b_maximum_per_year = .3 # [m/yr]

# Set this up to automatically number IceFlow outputs using glob
ic.output_filename=None
ic.output_figure=None
ic.plot_at_end_flag=False
ic.plot_during_run_flag = True
#ic.plot_t_years = ic.run_length_years
ic.boundary_condition = 'Dirichlet0'

ic.GRASS_raster_ice_extent = 'FakeMeasuredExtents'

ic.verbose = False

# Flexure
ic.isostatic = False

ic.initialize()
ic.run()
ic.finalize()


