import core as iceflow
ic = iceflow.IceFlow()

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
ic.dt_years = 2.

ic.elevation = 'elev'

ic.n = 792000
ic.s = 240000
ic.w = 108000
ic.e = 868000
ic.dx = 5000
ic.dy = 5000

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

ic.ELA = None
ic.dbdz = None
ic.bcap = None

# Set this up to automatically number IceFlow outputs using glob
ic.output_filename=None
ic.output_figure=None
ic.plot_t_years = ic.run_length_years

self = ic

# INITIALIZE
self.initialize_define_region()
self.initialize_GRASS()
self.initialize_elevation_x_y_and_ice_grids_from_GRASS()
if self.mass_balance_parameterization == 'TP_PDD':
  self.initialize_climate_from_GRASS()
elif self.mass_balance_parameterization == 'ELA':
  self.basic_mass_balance_with_GRASS()
self.initialize_compute_variables()
self.initialize_output_lists()
self.initialize_runtime()
self.enforce_boundary_conditions()
self.initialize_sparse_array()
ic.Pa /= (1000. * ic.secyr)

# gFlex

# Region
wf=-600000
nf=1200000
ef=1600000
sf=-200000
self.gr.run_command('g.region', n=nf, s=sf, e=ef, w=wf)

import gflex

flex = gflex.F2D()

flex.Quiet = False

flex.Method = 'FD'
flex.PlateSolutionType = 'vWC1994'
flex.Solver = 'direct'

flex.g = 9.8 # acceleration due to gravity
flex.E = 65E10 # Young's Modulus
flex.nu = 0.25 # Poisson's Ratio
flex.rho_m = 3300. # MantleDensity
flex.rho_fill = 0. # InfiillMaterialDensity

flex.dx = ic.dx
flex.dy = ic.dy
flex.BC_W = '0Displacement0Slope' # west boundary condition
flex.BC_E = '0Displacement0Slope' # east boundary condition
flex.BC_S = '0Displacement0Slope' # south boundary condition
flex.BC_N = '0Displacement0Slope' # north boundary condition

#flex.Te = 3700.
Tegrid = self.garray.array()
Tegrid.read('TeIcelandModified')
flex.Te = np.array(Tegrid)

qs_start = self.garray.array()
flex.qs = np.zeros(qs_start.shape) # Just to get it going

flex.initialize()

# Isostatic response
topo_initial = ic.Zb.copy()

dt_flexure = 100 # years
t_flexure_update = 0

def isostatic_response(equilibrium_deflection, topo_now, topo_initial, time_scale, dt):
  dz = (equilibrium_deflection - (topo_now - topo_initial)) * dt / time_scale
  return dz

flexure_time_scale = 1E-20 # instantaneous

dt_flexure = 100 # years
t_flexure_update = 0

# Return to computational region 
self.gr.run_command('g.region', n=self.n, s=self.s, w=self.w, e=self.e, ewres=self.dx, nsres=self.dy)

# RUN
for self.t_i in self.t_all:
  if self.t_i * ic.dt_years >= t_flexure_update:
    t_flexure_update += dt_flexure
    #grass.run_command('g.region', n=nf, s=sf, e=ef, w=wf)
    # ASSUMING SAME GRID CELL SIZE!
    flex.qs[(self.s-sf)/self.dx:-(nf-self.n)/self.dy, (self.w-wf)/self.dx:-(ef-self.e)/self.dx] = (ic.H - ic.H0) * 917. * 9.8
    flex.run()
    #grass.run_command('g.region', n=n, s=s, e=e, w=w)
    try:
      print np.min(dz), np.max(dz)
    except:
      pass
  ic.Zb = topo_initial + flex.w[(self.s-sf)/self.dx:-(nf-self.n)/self.dx, (self.w-wf)/self.dx:-(ef-self.e)/self.dx]
  # Still do isostatic calculation every time step: it is quick.
  self.build_sparse_array()
  self.solve_sparse_equation()
  if self.t_i % self.plot_t_years == 0:
    self.plot()


