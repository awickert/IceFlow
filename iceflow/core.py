import numpy as np
from matplotlib import pyplot as plt
import scipy
from scipy.io import loadmat
from scipy.sparse import linalg
from scipy.interpolate import RectBivariateSpline
import os
import sys
import warnings

class IceFlow(object):

  def __init__(self):
    self.add_default_variable_values()
    # Need to put defaults here, but can't stop user from running the model
    # if gFlex is not installed
    try:
      self.start_gflex_and_assign_default_variables()
    except:
      pass
    
  def initialize(self):
    self.initialize_define_region()
    if self.useGRASS:
      self.initialize_GRASS()
      self.initialize_elevation_x_y_and_ice_grids_from_GRASS()
      if self.mass_balance_parameterization == 'TP_PDD':
        self.initialize_climate_from_GRASS()
      elif self.mass_balance_parameterization == 'ELA':
        self.basic_mass_balance_with_GRASS()
    else:
      print "No mass balance methods without GRASS created yet!"
    if self.isostatic:
      self.initialize_isostatic()
    self.initialize_compute_variables()
    self.initialize_output_lists()
    self.enforce_boundary_conditions()
    self.initialize_sparse_array()
    if self.plot_during_run_flag:
      if self.isostatic:
        plt.figure(1, figsize=(12,8))
      else:
        plt.figure(1)
      plt.show(block = False)

  def update(self):
    if self.verbose:
      # display current time step in command window [a]
      print 'model year:', '%10.1f' %(self.t_years[self.t_i])
    if self.isostatic:
      self.update_isostatic_response()
    self.build_sparse_array()
    self.solve_sparse_equation()
    self.output()
    self.t_i += 1

  def output_time(self):
    """
    Checks if the time step is appropriate to write output, generate plots, 
    compare models against data, or do any as-of-yet-unknown additional 
    operations
    """
    if self.t_years[self.t_i] < \
      (self.record_timesteps_years[self.record_index] + self.dt/self.secyr/2.) \
      and self.t_years[self.t_i] >= \
      (self.record_timesteps_years[self.record_index] - self.dt/self.secyr/2.):
      self.record_model_parameters()
      outbool = True
    else:
      outbool = False
    return outbool
    
  def output(self):
    """
    At selected time steps t_i, record various model parameters, and, if so 
    desired, plot.
    """
    if self.output_time():
      if not self.quiet:
        # display current time step in command window [a]
        print 'Recording snapshot at model year:', \
              '%10.1f' %(self.t_years[self.t_i])
      if self.GRASS_raster_ice_extent:
        self.compare_with_ice_extent_GRASS_raster()
        if self.verbose:
          print "Model extent outside data region [%]:", \
                self.ModelOutsideData_FractOfIceAreaFromData[-1] * 100
          print "Data extent outside model region [%]:", \
                self.DataOutsideModel_FractOfIceAreaFromData[-1] * 100
      if self.record_model_parameters_flag:
        self.record_model_parameters()
      if self.plot_during_run_flag:
        self.plot_during_run()
    
  def run(self):
    for self.ts in self.t:
      self.update()

  def finalize(self):
    if not self.quiet:
      print "***"
      print "Model complete."
      if not self.quiet:
        print self.t_years[-1], "years elapsed"
      print "***"
    if self.output_filename:
      self.save_output()
    if self.OutNameGRASS:
      self.saveGRASS()
    if self.plot_at_end_flag:
      self.plot_at_end()
    if self.plot_during_run_flag:
      plt.figure(1)
      plt.show(block = True)
      
  #def resample_qs_to_Te(self, option):
  #  """
  #  Resample the load array to the (larger) Te array, which also might be at
  #  a different scale
  #  """
  #  if option == '

  def update_isostatic_response(self):
    if self.t_years[self.t_i] >= self.t_flexure_update:
      self.t_flexure_update += self.flexure_recalculation_frequency
      if type(self.flex.Te) == np.ndarray:
        # Old subgrid-extraction method
        # Requires the same dx
        # self.flex.qs[(self.south-self.flex.south)/self.flex.dy: \
        #             -(self.flex.north-self.north)/self.flex.dy, \
        #             (self.west-self.flex.west)/self.flex.dx: \
        #             -(self.flex.east-self.east)/self.flex.dx ] = \
        #                                 (self.H - self.H0) * 917. * 9.8
        # New interpolation method -- any dx, any grid
        dHinterp_fcn = RectBivariateSpline( self.x, self.y, 
                                            (self.H - self.H0).transpose() )
        dHinterp = dHinterp_fcn(self.flex.xcenter, self.flex.ycenter).transpose()
        # Set everything to 0 outside of ice model domain
        dHinterp_mask = (self.flex.Y <= self.north) * \
                        (self.flex.Y >= self.south) * \
                        (self.flex.X <= self.east)  * \
                        (self.flex.X >= self.west)
        dHinterp *= dHinterp_mask
        self.flex.qs = dHinterp * 917. * 9.8
      else:
        self.flex.qs = (self.H - self.H0) * 917. * 9.8
      self.flex.run()
      self.flex.finalize()
    self.compute_isostatic_response()

  def compute_isostatic_response(self):
    # exponential for small ts b/c of driving gradient getting smaller
    # Old subgrid-extraction method
    # Requires the same dx
    #equilibrium_deflection = self.flex.w[(self.south-self.flex.south)/self.flex.dy: \
    #                                     -(self.flex.north-self.north)/self.flex.dy, \
    #                                     (self.west-self.flex.west)/self.flex.dx: \
    #                                     -(self.flex.east-self.east)/self.flex.dx]
    # New interpolation method -- any dx, any grid
    # Probably downscaling here -- coarser grid allowed due to low-pass
    # filtering of flexure
    winterp_fcn = RectBivariateSpline( self.flex.xcenter, self.flex.ycenter, \
                                       self.flex.w.transpose() )
    equilibrium_deflection = winterp_fcn(self.x, self.y).transpose()
    # THIS ASSUMES NO EROSION!!!!!!!!!!!!!!
    # AS IN, APPROACHES EQUILIBRIUM AS ZB-ZB_INITIAL
    self.dz = (equilibrium_deflection - \
               (self.Zb - self.Zb_initial)) * self.dt_years / \
               self.isostatic_response_time_scale
    self.Zb += self.dz

  def add_default_variable_values(self):
    ############
    # Run time #
    ############
    self.run_length_years = None
    self.t_start_years = None
    self.dt_years   = 0.25
    self.t_i = 0 # counter
    
    #############
    # VERBOSITY #
    #############
    self.verbose = False
    self.debug = False # Not used yet
    self.quiet = False # If true, will set self.debug and self.verbose False
                       # during initialize() step.
  
    ###################
    # Basic constants #
    ###################
    self.secyr = 3600*24*365.24 # conversion factor 1 [s/a]
    self.g = 9.81 # gravitational acceleration [m/s2]
    self.rho = 917. # glacier density [kg/m3]
    self.R = 8.314 # ideal gas constant [J/mol/K]
    
    ########################################
    # Constants that we may what to change #
    ########################################
    self.nGlen = 3 # flow law exponent [unitless]
    # cold/warm ice creep activation energy [J/mol]
    self.QcCold = 60000.
    self.QcWarm = 139000.
    # cold/warm ice reference flow law parameter [/Pa3/a]
    self.AoCold = 1.14E-5
    self.AoWarm = 5.47E10
    
    ###################
    # Run with gFlex? #
    ###################
    self.isostatic = False
    # Using an elastic lithosphere with a response time scale
    # See paper showing that this is the best approximation to a true global model!!!!!!!!!!!!!!
    # <-----------------------------------------------------------------------------TO DO
    # Response time scale
    self.isostatic_response_time_scale = 4000. # [yr]
    # Elastic thickness -- string (filename or GRASS raster) or float
    self.ElasticThickness = None # [m]
    # How often to recalculate flexural response
    self.flexure_recalculation_frequency = 100 # [yr]
    # First time to apply flexural correction (this becomes a condition
    # for the update and changes through the run)
    self.t_flexure_update = 0 # [yr]
    
    #############
    # Plotting? #
    #############
    self.plot_during_run_flag = False # Plots at times when outputs are recorded
                                      # so depends on self.record_frequency_years
    self.plot_at_end_flag = True
    
    ##########
    # Output #
    ##########
    # Interval over which to record output [years]
    self.record_frequency_years = None
    # Interval over which to compare against data [years]
    # Will default to the same interval
    # And currentlly does not matter
    self.data_comparison_frequency_years = None # CURRENTLY NOT IMPLEMENTED!
                                                # can set it but data comparison
                                                # time interval tied to
                                                # self.record_frequency_years
                                                # (and this defaults to that
                                                # as well)
    # Flag to save model parameters
    self.record_model_parameters_flag = False
    # Flag to save misfits from data-model comparison
    self.data_comparison_frequency_years_flag = False
    # Output filename -- set to None if you do not want to save ouptput
    self.output_filename='IceFlow_Output'
    # Leave blank for no output figure
    self.output_figure = 'IceThickness.png'
    # GRASS GIS output?
    self.OutNameGRASS = None

    ###########################
    # Basal surface elevaiton #
    ###########################
    # If self.Zb_initial = a grid, then it is the bed elevaiton
    # If it is a string, then it is a file path or GRASS raster map
    # (and this is defined by useGRASS, bellow)
    self.elevation = None # MUST BE DEFINED IN INPUT!
    self.dx = None # Defined if run outside GIS
    self.dy = None # Will default to self.dx if it remains None and dx
                   # is defined
    # These next two are defined for the case in which there is no GRASS
    # GIS integration, to set the grid coordinate system
    self.xmin = None
    self.ymin = None
    
    ###############
    # GRASS setup #
    ###############
    self.useGRASS=False
    self.location=None
    self.mapset='PERMANENT'
    self.gisbase=None
    self.grassdata_dir='grassdata'
    
    #################
    # sliding setup #
    #################
    self.sliding=True
    # BASAL SLIDING COEFFICIENT [m/a/Pa]
    # Proportional to driving stress. Marshall et al. (2005) use 0.0006 m/a/Pa 
    # for Vatnajokull Ice Cap.
    # Colgan used 0.0003 m/a/Pa for the Front Range in Colorado.
    # Defaults to 0.0006 m/a/Pa
    self.C0_per_year=0.0006
    
    #############
    # Ice setup #
    #############
    # ICE TEMPERATURE [K]
    # Defaults to T=273.15 K: isothermal and constant through time and space.
    # Right now, is constant, and cannot be used with arrays
    # So flow law parameter defined here alone
    self.T=273.15 # Can also be a grid and/or vary with time
    # BOUNDARY CONDITIONS
    # The domain boundary should not intersect
    # any major ice masses. These default boundary condition options are suited
    # for dealing with small peripheral ice masses adjacent to the domain
    # boundary.
    # 
    # Dirichlet0 - prescribed head: boundary cells are prescribed as zero ice 
    # thickness. "Eliminated" flux is recorded as a dynamic leackage term.
    # 
    # Neumann0 - prescribed flux: boundary cells have zero flux leaving domain.
    # Ice "piles up" where it flows against the domain boundary.
    self.boundary_condition='Dirichlet0'
    self.StartingIce = None # Set to the name of a raster grid of initial ice 
                             # thicknesses,
                             # e.g., from a former model run or field data
    
    #################################
    # Mass balance setup -- generic #
    #################################
    # Maximum mass balance -- this can be a grid or a scalar
    self.b_maximum_per_year = 1.0
    # PICK PARAMETERIZATION TYPE
    # Options:
    #   TP_PDD --> Temperature and Precipitation, PDD
    #   ELA    --> Prescribed mass balance with ELA
    #              (constant if inputs are numbers, variable in space if inputs 
    #               are grids or scalars pointing to GRASS maps)
    self.mass_balance_parameterization = None
    
    ################################
    # Mass balance setup -- TP_PDD #
    ################################
    # Temperature and precipitation: can be modern or from the past
    # Strong to grid location, load value, etc.
    self.temperature_melt_season = None
    self.precipitation_solid     = None
    # If we have modern temperature and precipitation fields to modify
    # These could be scalar or grid values
    # And could (in theory) change with time, though they don't now.
    self.T_correction = 0. # air temperature offset to glacial climate [K]
    self.P_factor = 1. # precipitation scaling factor to glacial climate [fraction]
    # The following are used to project temperature and precipitation to altitude
    self.temperature_lapse = -4.7/1000. # See Anderson et al. (2014) for a good lapse rate compilation
    # precip_lapse doesn't take drying at high altitudes into account
    self.precip_lapse_years = 0.3698/1000. # characteristic regional precip. lapse rate [m/a/m]
    self.melt_factor_days = 6./1000. # surface ablation per air temperature (melt factor) [m/d/K]
    # Length of melt season
    self.melt_season_length_days = None

    #############################
    # Mass balance setup -- PMB #
    #############################
    # These can be grids or constant values across the whole system
    self.ELA  = None
    self.dbdz_per_year = None
    # (see above for b_maximum_per_year)

    ###################
    # Data comparison #
    ###################
    self.compare_at_intermediate_time_steps = True 
    self.GRASS_raster_ice_extent = None # Name of GRASS raster map with ice
                                        # extent from data
    # Very preliminary ideas on how to do the rest
    self.GRASS_vector_ice_extent = None # Not yet implemented
    self.GRASS_moraine_ages = None # Not yet implemented
    self.GRASS_glacial_polish_ages = None # Not yet implemented.
    

  def initialize_define_region(self):
    # Glacier model domain -- must be projected coordinate system
    # Should be defined after a grid is defined
    # These bounds work for both GRASS and non-GRASS initializations
    if self.north and self.south:
      pass
    else:
      if self.ymin is None:
        if not self.quiet:
          print "Automatically setting undefined south value to 0"
        self.ymin = 0
      self.south = self.ymin + self.dy/2.
      self.north = self.ymin + self.dy * self.Zb_initial.shape[0] - self.dy/2.
      self.y = np.arange(self.south, self.north + self.dy/10., self.dy)
    if self.west and self.east:
      pass
    else:
      if self.xmin is None:
        self.xmin = 0
        if not self.quiet:
          print "Automatically setting undefined west value to 0"
      self.west = self.xmin + self.dx/2.
      self.east = self.xmin + self.dx * self.Zb_initial.shape[1] - self.dy/2.
      self.x = np.arange(self.west, self.east + self.dx/10., self.dx)
    if self.west and self.east and self.north and self.south:
      pass
    else:
      self.X, self.Y = np.meshgrid(self.x, self.y) # And these are cell centers

  def initialize_compute_variables(self):
    # Variable conversions from years to seconds
    self.C0 = self.C0_per_year / self.secyr
    self.b_maximum = self.b_maximum_per_year / self.secyr
    self.dbdz = self.dbdz_per_year = self.secyr
    self.dt = self.dt_years*self.secyr # time step [s]
    # Combine cold and warm ice flow law parameters [/Pa3/a]
    self.A = self.AoCold*np.exp(-self.QcCold/self.R/self.T) * (self.T < 263.5) + self.AoWarm*np.exp(-self.QcWarm/self.R/self.T) * (self.T >= 263.5) / self.secyr
    # Time
    if self.record_frequency_years is None:
      #self.record_frequency_years = np.min((self.run_length_years, 10.))
      self.record_frequency_years = self.run_length_years / 10
      if self.record_frequency_years == 0:
        self.record_frequency_years = None # Then just don't record if it is 
                                           # a run this short -- or at least
                                           # make the user set the record 
                                           # frequency
    if self.data_comparison_frequency_years is None \
        and self.data_comparison_frequency_years is not None:
      # I guess the "not none" is redundant
      self.data_comparison_frequency_years = self.record_frequency_years

    else:
      pass # must have been defined by user
    tEnd = self.run_length_years*self.secyr # model run length [s]
    self.t = np.arange(0, tEnd+self.dt/10., self.dt) # time vector [s]
    self.t_years = self.t / self.secyr
    self.record_timesteps_years = np.arange(0, self.run_length_years+self.record_frequency_years/2., self.record_frequency_years)
    # Input elevation map
    if type(self.elevation) is not str and type(self.elevation):
      self.Zb_initial = self.elevation
    if self.Zb_initial is not None:
      self.Zb = self.Zb_initial.copy()
    # Second thoughts on this
    #if self.dx and self.dy is None:
    #  self.dy = self.dx.copy()
    # PDD
    # Surface mass balance - An elevation-dependent surface mass 
    # balance is parameterized from observed contemporary solid precipitation 
    # and air temperature fields that are perturbed to glacial climate.
    self.precip_lapse = self.precip_lapse_years/self.secyr # convert units of precipitation lapse rate [m/s/m]
    self.mu = self.melt_factor_days/86400. # convert units of melt factor [m/s/K]
    self.melt_season_length = self.melt_season_length_days * 86400.
    # Verbosity
    if self.quiet:
      self.verbose = False
      self.debug = False
    # output data--model comparison lists
    if self.GRASS_raster_ice_extent:
      self.ModelOutsideData_FractOfIceAreaFromData = []
      self.DataOutsideModel_FractOfIceAreaFromData = []
    # Plotting with gFlex
    self.dz = np.zeros(self.Zb.shape) # define dz for first plotting step
    # grid
    try:
      self.x
    except:
      self.x = self.X[0,:]
    try:
      self.y
    except:
      self.y = np.flipud(self.Y[:,0])
    
  def initialize_output_lists(self):
    self.record_index = 0 # index of recorded selected time steps [unitless]
    self.b_timestep = 0 # surface mass balance record counter [a]
    self.a_timestep = 0 # surface ablation record counter [a]
    self.c_timestep = 0 # surface accumulation record counter [a]
    self.outflow_timestep = 0 # domain outflow record counter [a]
    self.time_series = np.zeros((len(self.record_timesteps_years),5)) # record of selected variables at selected time step [variable]  
    self.H_record = [] # List to record ice thicknes at time slices [m]
    self.Zb_record = [] # List to record bed elevation at time slices [m] -- important for isostasy, maybe in future if erosion is included
    self.dz_record = [] # change in bed elev (isostasy, etc.)
    self.uS_record = [] # List to record sliding velocities at time slices [m/a]
    self.b_record = [] # List to record mass balance at time slices [m/a]
    self.uD_record = [] # List to record depth-averaged deformational velocities [m/a]
    self.t_record = [] # time step [yr]

  def initialize_GRASS(self):
    if self.location:
      if not self.quiet:
        print "Using GRASS location", self.location
    else:
      sys.exit("Must define a GRASS GIS location.")
    # Implicitly 
    os.environ['GISBASE'] = self.gisbase
    gisdbase = os.path.join(os.environ['HOME'], self.grassdata_dir)
    sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
    import grass.script as gr
    import grass.script.setup as gsetup
    gsetup.init(self.gisbase, gisdbase, self.location, self.mapset)
    from grass.script import array as garray
    self.gr = gr
    self.garray = garray
    if self.north is None and self.south is None and self.west is None and self.east is None:
      pass
    elif self.north is None or self.south is None or self.west is None or self.east is None:
      sys.exit('Incomplete description of edge values')
    else:
      self.gr.run_command('g.region', n=self.north, s=self.south, w=self.west, e=self.east)
      # Then apply the desired resolution.
      # This won't be *quite* a square grid, but will be close.
      if self.dx and self.dy:
        self.gr.run_command('g.region', ewres=self.dx, nsres=self.dy)
      else:
        if not self.quiet:
          print "dx and dy not both defined: not updating grid resolution."
    self.gr.run_command('g.region', save='IceFlowRegion', overwrite=True) # Save it
    
  def initialize_elevation_x_y_and_ice_grids_from_GRASS(self):
    """
    Get elevation grids from GRASS as well as the grid dimensions
    and the pre-existing ice (if applicable)
    """
    # Import arrays of topography
    self.Zb_initial = self.garray.array()
    self.Zb_initial.read(self.elevation)
    self.Zb = self.Zb_initial.copy() # in case isostasy is used
    # And then create the grids of the x and y values that go along with these positions
    self.dx = self.gr.region()['ewres']
    self.dy = self.gr.region()['nsres']
    self.gr.mapcalc('x = x()', overwrite=True, quiet=True) # is there a better way?
    self.gr.mapcalc('y = y()', overwrite=True, quiet=True) # is there a better way?
    #X = np.arange(0, Z.shape[1], self.dx)
    #Y = np.arange(0, Z.shape[0], self.dy)
    self.X = self.garray.array() # cell centers
    self.Y = self.garray.array() # cell centers
    self.X.read('x') # cell centers
    self.Y.read('y') # cell centers
    self.ny, self.nx = self.X.shape # number of easting and northing nodes [unitless]  
    self.H0 = self.garray.array() # Starts out with 0's at the right shape
    # And ice array, if applicable
    if self.StartingIce:
      self.H0.read(self.StartingIce)
    else:
      self.H0[...] = 0
    # At start, H = H0
    self.H = self.H0.copy()

  def initialize_climate_from_GRASS(self):
    """
    Get grids from GRASS
    Run after elevation grids obtained (this needs to be fixed)
    """
    self.Tair = self.garray.array()
    self.Tair.read(self.temperature_melt_season)
    self.Pair = self.garray.array()
    self.Pair.read(self.precipitation_solid)
    # NOTE:
    # self.Pair /= (1000. * self.secyr) # e.g., PRISM is in mm/yr, change to m/s.
    # Correction factors
    self.Ta = self.Tair + self.T_correction # Temperature average through melt season [degC]
    self.Pa = self.Pair * self.P_factor # Total annual precipitation (assume snow?) (doesn't have to be this) [m/s]

  def basic_mass_balance_with_GRASS(self, ELA, dbdz):
    """
    dbdx = change in mass balance w/ change in x
    dbdy = ditto for y
    ...
    """
    self.ela0 = self.ELA.copy() # Initial ELA field -- may change over time

  def enforce_boundary_conditions(self):
    # boundary conditions:  
    if self.boundary_condition == 'Dirichlet0':
      # prescribed zero ice thickness
      # boundary condition mask [binary]
      self.BC = np.ones((self.ny,self.nx))
      self.BC[:1,:] = 0 
      self.BC[:,:1] = 0 
      self.BC[-1:,:] = 0 
      self.BC[:,-1:] = 0
      self.Zb_initial *= self.BC
    elif self.boundary_condition == 'Neumann0':
      # type 2 (prescribed zero ice flux) 
      self.BC = np.ones((self.ny,self.nx)) # boundary condition mask [binary]
    
  def initialize_sparse_array(self):

    self.R_term_yes = np.hstack((np.ones((self.ny,self.nx-1)), np.zeros((self.ny,1)))) # identify nodes where the RIGHT-hand matrix term is present [binary]
    self.L_term_yes = np.hstack((np.zeros((self.ny,1)), np.ones((self.ny,self.nx-1)))) # identify nodes where the LEFT-hand matrix term is present [binary]
    self.D_term_yes = np.vstack((np.ones((self.ny-1,self.nx)), np.zeros((1,self.nx)))) # identify nodes where the DOWN-ward matrix term is present [binary]  
    self.U_term_yes = np.vstack((np.zeros((1,self.nx)), np.ones((self.ny-1,self.nx)))) # identify nodes where the UP-ward matrix term is present [binary]  
    self.C_term_yes = np.ones((self.ny,self.nx)) # identify nodes where the CENTRE matrix term is present [binary]
       
    # bedrock/land surface field:
    #self.Zb = interp2(XI_500m,YI_500m,Z_500m,x,y);
    #  # resample bedrock elevation at nodes j,i [m] -- not doing here b/c should happen before this code -- Landlab integration
    self.Zb *= self.BC # constrain bedrock elevation [m]
    #
    # MUST MOVE if changing elevation (e.g., flexure, erosion, asteroid impacts, global thermonuclear war...)
    self.ZbP1 = np.hstack((self.Zb[:,1:], np.zeros((self.ny,1)))) # bedrock elevation at node j,i+1 [m]
    self.ZbM1 = np.hstack((np.zeros((self.ny,1)), self.Zb[:,:-1])) # bedrock elevation at node j,i-1 [m]
    self.Zb_jP1i = np.vstack((self.Zb[1:,:], np.zeros((1,self.nx)))) # bedrock elevation at node j+1,i [m]
    self.Zb_jM1i = np.vstack((np.zeros((1,self.nx)), self.Zb[:-1,:])) # bedrock elevation at node j-1,i [m]  
      
    # initialize evolving ice geometry variables:
    self.H = self.H0.copy() # ice thickness at node j,i [m] 
    self.Zs = self.Zb + self.H # ice surface elevation at node j,i [m]
    
  def build_sparse_array(self):
    # calculate ice form and flow variables:
    self.Zs = self.Zb + self.H # ice surface elevation at node j,i [m]
    if self.mass_balance_parameterization == 'TP_PDD':
      self.dz_plus_H = self.H + self.Zb - self.Zb_initial
      a = (self.Ta + (self.dz_plus_H*self.temperature_lapse))*self.mu * self.melt_season_length # surface ablation scaled to melt season [m/yr]
      c = (self.Pa + (self.dz_plus_H*self.precip_lapse) )*self.secyr # surface accumulation [m/yr]
      self.b = (c - a) # surface mass balance [m/s]
    elif self.mass_balance_parameterization == 'ELA':
      self.b = self.dbdz*(self.Zs-self.ela0)
    self.b[self.b > self.b_maximum_per_year] = self.b_maximum_per_year
    self.b /= self.secyr # No clue whether I was trying to do this or not -- looks like seconds are a go!
      
    self.HPh = np.hstack(( (self.H[:,1:]+self.H[:,:-1])/2., np.zeros((self.ny,1)) )) # ice thickness at node j,i+1(?) [m]
    self.HMh = np.hstack(( np.zeros((self.ny,1)), (self.H[:,1:] + self.H[:,:-1])/2. )) # ice thickness at node j,i-1(?) [m]
    self.H_jPhi = np.vstack(( (self.H[1:,:]+self.H[:-1,:])/2., np.zeros((1,self.nx)) )) # ice thickness at node j+1(?),i [m]
    self.H_jMhi = np.vstack(( np.zeros((1,self.nx)), (self.H[1:,:]+self.H[:-1,:])/2. )) # ice thickness at node j-1(?),i [m]      
    
    self.alpha = (np.hstack(( np.zeros((self.ny,1)), (self.Zs[:,2:] - self.Zs[:,:-2]) / (self.X[:,2:] - self.X[:,:-2]), np.zeros((self.ny,1)) ))**2 + np.vstack(( np.zeros((1,self.nx)), (self.Zs[2:,:] - self.Zs[:-2,:])/(self.Y[2:,:] - self.Y[:-2,:]), np.zeros((1,self.nx)) ))**2)**0.5 # absolute ice surface slope at node i,j [m/m]
    self.alphaPh = np.hstack(( (self.alpha[:,1:]+self.alpha[:,:-1])/2., np.zeros((self.ny,1)) )) # absolute ice surface slope at node j,i+1(?) [m]
    self.alphaMh = np.hstack(( np.zeros((self.ny,1)), (self.alpha[:,1:]+self.alpha[:,:-1])/2. )) # absolute ice surface slope at node j,i-1(?) [m]
    self.alpha_jPhi = np.vstack(( (self.alpha[1:,:]+self.alpha[:-1,:])/2., np.zeros((1,self.nx)) )) # absolute ice surface slope at node j+1(?),i [m]
    self.alpha_jMhi = np.vstack(( np.zeros((1,self.nx)), (self.alpha[1:,:]+self.alpha[:-1,:])/2. )) # absolute ice surface slope at node j-1(?),i [m]      
        
    self.dZsdxPh = np.hstack(( (self.Zs[:,1:]-self.Zs[:,:-1]) / (self.X[:,1:]-self.X[:,:-1]), np.zeros((self.ny,1)) )) # directional ice surface slope at node j,i+1(?) [m/m]
    self.dZsdxMh = np.hstack(( np.zeros((self.ny,1)), (self.Zs[:,1:]-self.Zs[:,:-1]) / (self.X[:,1:]-self.X[:,:-1]) )) # directional ice surface slope at node j,i-1(?) [m/m]
    self.dZsdy_jPhi = np.vstack(( (self.Zs[1:,:]-self.Zs[:-1,:]) / (self.Y[1:,:]-self.Y[:self.ny-1,:]), np.zeros((1,self.nx)) )) # directional ice surface slope at node j+1(?),i [m/m]
    self.dZsdy_jMhi = np.vstack(( np.zeros((1,self.nx)), (self.Zs[1:,:]-self.Zs[:-1,:]) / (self.Y[1:,:]-self.Y[:-1,:]) )) # directional ice surface slope at node j-1(?),i [m/m]
      
    tauPh = -self.rho*self.g*self.HPh*self.dZsdxPh # driving stress at node j,i+1(?) [Pa] - positive x direction
    tauMh = -self.rho*self.g*self.HMh*self.dZsdxMh # driving stress at node j,i-1(?) [Pa] - negative x direction
    tau_jPhi = -self.rho*self.g*self.H_jPhi*self.dZsdy_jPhi # driving stress at node j+1(?),i [Pa] - positive y direction
    tau_jMhi = -self.rho*self.g*self.H_jMhi*self.dZsdy_jMhi # driving stress at node j-1(?),i [Pa] - negative y direction
    
    qPh = 2*self.A/(self.nGlen+2)*(self.rho*self.g*self.alphaPh)**(self.nGlen-1)*self.HPh**(self.nGlen+1)*tauPh # ice discharge at node j,i+1(?) [m2/s]
    qMh = 2*self.A/(self.nGlen+2)*(self.rho*self.g*self.alphaMh)**(self.nGlen-1)*self.HMh**(self.nGlen+1)*tauMh # ice discharge at node j,i-1(?) [m2/s]
    q_jPhi = 2*self.A/(self.nGlen+2)*(self.rho*self.g*self.alpha_jPhi)**(self.nGlen-1)*self.H_jPhi**(self.nGlen+1)*tau_jPhi # ice discharge at node j+1(?),i [m2/s]
    q_jMhi = 2*self.A/(self.nGlen+2)*(self.rho*self.g*self.alpha_jMhi)**(self.nGlen-1)*self.H_jMhi**(self.nGlen+1)*tau_jMhi # ice discharge at node j-1(?),i [m2/s]
    
    self.uD = ( ((qPh/np.maximum(self.HPh,1E-8) + qMh/np.maximum(self.HMh,1E-8)) / 2.)**2 + ((q_jPhi/np.maximum(self.H_jPhi,1E-8) + q_jMhi/np.maximum(self.H_jMhi,1E-8)) / 2.)**2)**0.5 # absolute depth-averaged deformational velocity at node j,i [m/s]

    if self.sliding:
      uSPh = -self.C0*tauPh # basal sliding velocity at node j,i+1(?) [m/s]
      uSMh = -self.C0*tauMh; 
        # basal sliding velocity at node j,i-1(?) [m/s]
      uS_jPhi = self.C0*tau_jPhi # basal sliding velocity at node j+1(?),i [m/s] # (AW) made these positive, solution looks better now
      uS_jMhi = self.C0*tau_jMhi # basal sliding velocity at node j-1(?),i [m/s] # (AW) but haven't checked why
      self.uS = (((uSPh + uSMh)/2)**2 + ((uS_jPhi + uS_jMhi)/2)**2)**0.5 # absolute basal sliding velocity at node j,i [m/s]       
    else:
      uSPh = uSMh = uS_jPhi = uS_jMhi = self.uS = 0
      
    self.dZsdxPh[self.dZsdxPh == 0] = 1E-6 # directional ice surface slope at node j,i+1(?) [m/m]
    self.dZsdxMh[self.dZsdxMh == 0] = 1E-6 # directional ice surface slope at node j,i-1(?) [m/m]
    self.dZsdy_jPhi[self.dZsdy_jPhi == 0] = 1E-6 # directional ice surface slope at node j+1(?),i [m/m]
    self.dZsdy_jMhi[self.dZsdy_jMhi == 0] = 1E-6 # directional ice surface slope at node j-1(?),i [m/m]
      
    DPh = qPh/self.dZsdxPh; DPh[:,-1] = 0 # diffusion term at node j,i+1(?) [m2/s]
    DMh = qMh/self.dZsdxMh; DMh[:,0] = 0 # diffusion term at node j,i-1(?) [m2/s]     
    D_jPhi = q_jPhi/self.dZsdy_jPhi; D_jPhi[-1,:] = 0 # diffusion term at node j+1(?),i [m2/s]
    D_jMhi = q_jMhi/self.dZsdy_jMhi; D_jMhi[0,:] = 0 # diffusion term at node j-1(?),i [m2/s]  
    
    # create and solve 5-banded matrix:  
    Array_L = (+ DMh*self.dt/self.dx**2 + uSMh*self.dt/2/self.dx) * self.L_term_yes # unknown at node j,i-1 [unitless]
    Array_D = (+ D_jMhi*self.dt/self.dy**2 + uS_jMhi*self.dt/2/self.dy) * self.D_term_yes # unknown at node j-1,i [unitless]
    Array_C = 1 \
      + (- DPh*self.dt/self.dx**2 - uSPh*self.dt/2./self.dx) * self.R_term_yes \
      + (- DMh*self.dt/self.dx**2 + uSMh*self.dt/2./self.dx) * self.L_term_yes \
      + (- D_jPhi*self.dt/self.dy**2 - uS_jPhi*self.dt/2./self.dy) * self.U_term_yes \
      + (- D_jMhi*self.dt/self.dy**2 + uS_jMhi*self.dt/2./self.dy) * self.D_term_yes
      # unknown at node j,i [unitless]
    Array_U = (+ D_jPhi*self.dt/self.dy**2 - uS_jPhi*self.dt/2/self.dy) * self.U_term_yes # unknown at node j+1,i [unitless]
    Array_R = (+ DPh*self.dt/self.dx**2 - uSPh*self.dt/2/self.dx) * self.R_term_yes # unknown at node j,i+1 [unitless]
    Array_rhs = (+ self.b*self.dt + self.H) * self.C_term_yes \
      + DPh*(self.Zb - self.ZbP1)*self.dt/self.dx**2 * self.R_term_yes \
      - DMh*(self.ZbM1 - self.Zb)*self.dt/self.dx**2 * self.L_term_yes \
      + D_jPhi*(self.Zb - self.Zb_jP1i)*self.dt/self.dy**2 * self.U_term_yes \
      - D_jMhi*(self.Zb_jM1i - self.Zb)*self.dt/self.dy**2 * self.D_term_yes # known (right hand side vector) at node j,i [m] 
      
    Vec_D = np.reshape(Array_D, -1, order='F') # reshape unknown at node j-1,i from array to vector [unitless]     
    Vec_L = np.reshape(Array_L, -1, order='F') # reshape unknown at node j,i-1 from array to vector [unitless]   
    Vec_C = np.reshape(Array_C, -1, order='F') # reshape unknown at node j,i from array to vector [unitless]   
    Vec_R = np.reshape(Array_R, -1, order='F') # reshape unknown at node j,i+1 from array to vector [unitless]   
    Vec_U = np.reshape(Array_U, -1, order='F') # reshape unknown at node j+1,1 from array to vector [unitless]   

    # NECESSARY if updating these parameters
    self.ZbP1 = np.hstack((self.Zb[:,1:], np.zeros((self.ny,1)))) # bedrock elevation at node j,i+1 [m]
    self.ZbM1 = np.hstack((np.zeros((self.ny,1)), self.Zb[:,:-1])) # bedrock elevation at node j,i-1 [m]
    self.Zb_jP1i = np.vstack((self.Zb[1:,:], np.zeros((1,self.nx)))) # bedrock elevation at node j+1,i [m]
    self.Zb_jM1i = np.vstack((np.zeros((1,self.nx)), self.Zb[:-1,:])) # bedrock elevation at node j-1,i [m]  

    self.Vec_rhs = np.reshape(Array_rhs, -1, order='F') # reshape known at node j,i from array to vector [unitless]   
    
    # Sparse matrix
    diags = np.vstack(( np.hstack(( Vec_L[self.ny:], np.zeros(self.ny) )), \
                        np.hstack(( Vec_D[1:], 0 )), \
                        Vec_C, \
                        np.hstack(( 0, Vec_U[:-1] )), \
                        np.hstack(( np.zeros(self.ny), Vec_R[:-self.ny] )) ))
    self.Matrix_lhs = scipy.sparse.spdiags(diags, [-self.ny,-1,0,1,self.ny], self.ny*self.nx, self.ny*self.nx, format='csr') # create five-banded sparse matrix [unitless]

  def solve_sparse_equation(self):
    # solve ice thickness, in Matlab using backslash, here using umfpack [m]
    Vec_H = linalg.spsolve(self.Matrix_lhs, self.Vec_rhs, use_umfpack=False)

    self.H = np.reshape(Vec_H, (self.nx, self.ny)).transpose() # reshape ice thickness from vector to array [m]
    self.H[self.H < 1E-8] = 0 # constrain potentially negative ice thicknesses [m]
    
    # potentially implement domain boundary condition:
    H_pre = self.H # note ice thickness before implementing boundary condition [m]
    self.H = self.H*self.BC # implement prescribed boundary condition on ice thickness [m]
    H_post = self.H # note ice thickness again after implementing boundary condition [m]
    
  def update_recorded_mass_balance(self):
    # update mass balance terms between recorded time steps:
    self.outflow_timestep = self.outflow_timestep + sum(sum(H_pre - H_post))* \
      self.dx*self.dy*self.rho # update domain outflow (kg/a)  
    self.b_timestep = self.b_timestep + np.sum(self.b*(self.H>0))*self.dx*self.dy*self.dt*self.rho # update time step surface mass balance (kg/a)
    self.a_timestep = self.a_timestep + np.sum(a*(self.H>0))*self.dx*self.dy*self.dt*self.rho # update time step surface ablation (kg/a)
    self.c_timestep = self.c_timestep + np.sum(c*(self.H>0))*self.dx*self.dy*self.dt*self.rho # update time step surface accumulation (kg/a) 
  
  def compare_with_ice_extent_GRASS_raster(self):
    iceExtentData = self.garray.array()
    iceExtentData.read(self.GRASS_raster_ice_extent) # Must be binary mask of 1 (ice) and 
                                          # 0 (not ice)
    modelice = self.H > 1. # Just set a 1-meter limit to the ice thickness
                           # to define ice extent
                           # Might use slope-based subgrid methods in the future
    comparison = modelice + iceExtentData
    non_overlap = np.sum(comparison == 1) # number of cells that are ice in the 
                                          # model but are not in the data or
                                          # vice versa (are ice in the data 
                                          # but not in the model
    overlap_grid = modelice * iceExtentData
    model_farther_than_ice_margins = np.sum(modelice - overlap_grid)
    ice_margins_farther_than_model = np.sum(iceExtentData - overlap_grid)
    measured_ice_extent_cells = np.sum(iceExtentData) # could take this outside
    self.ModelOutsideData_FractOfIceAreaFromData.append( \
         model_farther_than_ice_margins/float(measured_ice_extent_cells))
    self.DataOutsideModel_FractOfIceAreaFromData.append( \
         ice_margins_farther_than_model/float(measured_ice_extent_cells))
  
  def start_gflex_and_assign_default_variables(self):
    import gflex
    self.flex = gflex.F2D()
    self.flex.Quiet = False
    self.flex.Method = 'FD'
    self.flex.PlateSolutionType = 'vWC1994'
    self.flex.Solver = 'direct'
    self.flex.g = 9.8 # acceleration due to gravity
    self.flex.E = 65E9 # Young's Modulus
    self.flex.nu = 0.25 # Poisson's Ratio
    self.flex.rho_m = 3300. # MantleDensity
    self.flex.rho_fill = 0. # InfiillMaterialDensity
    self.flex.BC_W = '0Displacement0Slope' # west boundary condition
    self.flex.BC_E = '0Displacement0Slope' # east boundary condition
    self.flex.BC_S = '0Displacement0Slope' # south boundary condition
    self.flex.BC_N = '0Displacement0Slope' # north boundary condition
    self.flex.north = None
    self.flex.south = None
    self.flex.west = None
    self.flex.east = None
    self.flex.dx = None
    self.flex.dy = None

  def initialize_isostatic(self):
    """
    Get gFlex going with standard selections for variables
    (can be overwritten by user)
    """
    if type(self.ElasticThickness) == str:
      if self.useGRASS:
        try:
          self.gr.run_command('g.region', rast=self.ElasticThickness)
          reg = self.gr.region()
          if self.flex.dx is None:
            self.flex.dx = reg['ewres']
          if self.flex.dy is None:
            self.flex.dy = reg['nsres']
          self.gr.run_command('g.region', n=self.flex.north, s=self.flex.south, \
                              w=self.flex.west, e=self.flex.east, \
                              nsres=self.flex.dy, ewres=self.flex.dx )
          self.gr.run_command('g.region', save='IceFlow_gFlexTe', overwrite=True)
          Tegrid = self.garray.array()
          Tegrid.read(self.ElasticThickness)
          self.flex.Te = np.array(Tegrid)
          self.gr.run_command('g.region', region='IceFlowRegion')
        except:
          pass
      else:
        sys.exit("Variable Te only works with GRASS")
    elif (type(self.ElasticThickness) == float) or \
         (type(self.ElasticThickness) == int):
      self.flex.Te = self.ElasticThickness
    else:
      sys.exit("No defined Te")
    if type(self.flex.Te) != np.ndarray:
      # For constant Te -- of they are not defined
      if self.flex.dx is None:
        self.flex.dx = self.dx
        self.flex.dy = self.dy
    
    # Starting loads -- assume 0
    if type(self.flex.Te) == np.ndarray:
      self.flex.qs = 0 * self.flex.Te
    else:
      self.flex.qs = 0 * self.H0
      
    # Grid on-cell-centers for Te
    self.flex.xcenter = np.arange(self.flex.west + self.flex.dx/2., \
                                  self.flex.east, self.flex.dx)
    self.flex.ycenter = np.arange(self.flex.south + self.flex.dy/2., \
                                  self.flex.north, self.flex.dy)
    self.flex.X, self.flex.Y = np.meshgrid(self.flex.xcenter, self.flex.ycenter)
    
    # Initialize at last
    self.flex.initialize()

  def record_model_parameters(self):
    self.record_timesteps_years = np.arange(0, self.run_length_years+self.record_frequency_years/2., self.record_frequency_years) # time-steps to be record [a]
    self.H_record.append(self.H)
      # record time step ice thickness field [m]
    self.Zb_record.append(self.Zb.copy())
      # record bed elevation [m/a] -- important if isostasy is implemented.
    self.dz_record.append(self.Zb - self.Zb_initial)
      # record bed elevation change from isostasy [m]
    self.uS_record.append(self.uS*self.secyr)
      # record time step basal sliding velocity field [m/a]
    self.uD_record.append(self.uD*self.secyr)
      # record time step deformational velocity field [m/a]
    self.b_record.append(self.b*self.secyr)
      # record mass balance [m/a]
    self.t_record.append(self.t[self.t_i]/self.secyr)
    #self.time_series[self.record_index,:] = np.hstack(( self.record_timesteps_years[self.record_index], \
    #                                          self.c_timestep/self.record_frequency_years, \
    #                                          self.a_timestep/self.record_frequency_years, \
    #                                          self.b_timestep/self.record_frequency_years, \
    #                                          self.outflow_timestep/self.record_frequency_years ))
      # update time series of mass balance elements [kg/a]
    self.record_index += 1;
      # update record index [unitless]
    self.b_timestep = 0;
      # reset surface mass balance counter [a]
    self.a_timestep = 0;
      # reset surface ablation counter [a]
    self.c_timestep = 0;
      # reset surface accumulation counter [a]
    self.outflow_timestep = 0;
      # reset domain outflow counter [a]

  def save_output(self):
    """
    save simulation output into a single .npy file that can be called on
    for graphical output at a later time
    
    Also save outputs of data--model comparison, if available
    """
    # State variables
    out_arrays = (self.H_record, self.Zb_record, self.uD_record, \
                  self.uS_record, self.b_record, self.t_record, \
                  self.time_series, self.T, self.A, self.C0, self.X, self.Y, \
                  self.Zb, self.BC, self.T_correction, \
                  self.P_factor, self.mu, self.dx, self.dy)
    np.save(output_filename, out_arrays)
    # Data -- model comparison
    data_model_comparison = (self.ModelOutsideData_FractOfIceAreaFromData, \
                             self.DataOutsideModel_FractOfIceAreaFromData )
    np.save(output_filename + '_DataComparison', data_model_comparison)
    
  def saveGRASS(self):
    iceGA = garray.array()
    iceGA[...] = self.H
    iceGA.write(self.OutNameGRASS+'_Hice', overwrite=True)
    iceGA[...] = self.uD
    iceGA.write(self.OutNameGRASS+'_uDeformation', overwrite=True)
    iceGA[...] = self.uS
    iceGA.write(self.OutNameGRASS+'_uSliding', overwrite=True)
    iceGA[...] = self.Zb
    iceGA.write(self.OutNameGRASS+'_zBed', overwrite=True)
    iceGA[...] = self.Zs
    iceGA.write(self.OutNameGRASS+'_zSurface', overwrite=True)

  def plot_at_end(self):
    """
    save == True: save figure images
    save == False: draw the plots on screen
    """

    plt.figure(2)
    plt.imshow(self.H, interpolation='nearest')
    # Catch warning about colorbar with all 0's in array
    with warnings.catch_warnings(record=True) as warning_messages:
      plt.colorbar()
    plt.title('Ice Thickness', fontsize=16)
    if self.output_figure:
      plt.savefig(self.output_figure)
      plt.close()
    plt.show()

  def plot_during_run(self):
    plt.clf()
    if self.isostatic:
      plt.subplot(221)
      plt.title('Ice thickness [m]')
      plt.imshow(self.H, interpolation='nearest')
      with warnings.catch_warnings(record=True) as warning_messages:
        plt.colorbar()
      plt.subplot(222)
      plt.title('Isostatic adjustment rate [m/yr]')
      plt.imshow(self.dz/self.dt_years, interpolation='nearest')
      with warnings.catch_warnings(record=True) as warning_messages:
        plt.colorbar()
      plt.subplot(223)
      plt.title('Surface mass balance [mm/yr]')
      plt.imshow(self.b * self.secyr * 1000., interpolation='nearest')
      with warnings.catch_warnings(record=True) as warning_messages:
        plt.colorbar()
      plt.subplot(224)
      plt.title('Total isostatic adjustment [m]')
      plt.imshow(self.Zb - self.Zb_initial, interpolation='nearest')
      with warnings.catch_warnings(record=True) as warning_messages:
        plt.colorbar()
    else:
      plt.imshow(self.H, interpolation='nearest')
      with warnings.catch_warnings(record=True) as warning_messages:
        plt.colorbar()
      plt.title('Ice Thickness', fontsize=16)
    plt.tight_layout()
    plt.draw()

