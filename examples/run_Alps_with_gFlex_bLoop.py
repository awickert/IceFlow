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

ic.run_length_years = 100. # years
ic.t_start_years = 0. # years
ic.dt_years = 5
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
ic.plot_during_run_flag = False
#ic.plot_t_years = ic.run_length_years
ic.boundary_condition = 'Dirichlet0'

ic.GRASS_raster_ice_extent = 'IceExtentAlpenLGM'
#ic.GRASS_raster_ice_extent = 'LGM_Alpen_noholes'

ic.verbose = False

ic.record_frequency_years = 20.

# Flexure
ic.isostatic = False
#ic.ElasticThickness = 20000
ic.ElasticThickness = 'EETEurope_meters' # Cell size is 20000 m; this is
                                         # minimum grid size for flexure
# Bounding box for flexure
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

# Some starting values just to get the model going
ic.ELA = 2300
ic.dbdz_per_year = 1E-3 # 1 m/yr over 1 km -- EDIT THIS TO FIT DATA!
ic.b_maximum_per_year = 0.3 # [m/yr]

# Initialize it
#ic.quiet=True
ic.initialize()
#ic.quiet=False

# Now make some lists of the variables that you want to change
# Get an empty array of the proper size
empty_array = np.array(ic.garray.array()) + 1.

# Then use this empty array with some rules to generate the gradients that
# you want to test.

# These look boring here, but can be fully spatially-distributed

ELAgrid1 = 2300 * empty_array
ELAgrid2 = 1500 * empty_array

BMAXgrid1 = 0.3 * empty_array
BMAXgrid2 = 0.8 * empty_array

# You can also use GIS plots to generate your inputs:
#BMAXgridFromGIS_1 = ic.garray.array()
#BMAXgridFromGIS_1.read('NameOfMyInputArray')
# And if you like, do mathematical operations on it
#BMAXgridFromGIS_1 *= .5

# Fill these brackets with arrays that you create
ELAgrids = [ELAgrid2]
BMAXgrids = [BMAXgrid2]

# And define lists of output values
ModelOutsideData_FractOfIceAreaFromData_list = []
DataOutsideModel_FractOfIceAreaFromData_list = []
RecordTimes_list = []

# Then run the code
for ic.ELA in ELAgrids:
  for ic.b_maximum_per_year in BMAXgrids:
    ic.reinitialize_mass_balance()
    print np.max(ic.Zb_initial)
    ic.run()
    ic.finalize()
    print np.max(ic.Zb_initial)
    # These output grids contain fractional areas of modeled ice extent that
    # lies outside the measured ice extent region and of measured ice extent
    # that lies outside the modeled ice extent region. Ideally both would be 0.
    ModelOutsideData_FractOfIceAreaFromData_list.append( \
      ic.ModelOutsideData_FractOfIceAreaFromData)
    DataOutsideModel_FractOfIceAreaFromData_list.append( \
      ic.DataOutsideModel_FractOfIceAreaFromData)
    # Times at which these misfits are recorded
    # With adaptive time stepping, these values may be off by a few months
    # to years, but this is a small issue right now
    RecordTimes_list.append(ic.record_timesteps_years)

# Consider adding something to save the grids with amount of overlap here
out_lists = (ModelOutsideData_FractOfIceAreaFromData_list, \
             DataOutsideModel_FractOfIceAreaFromData_list, \
             RecordTimes_list)
np.save('TestOut', out_lists) # Change the filename or make it a variable

