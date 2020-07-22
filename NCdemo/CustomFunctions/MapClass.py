""" 
A class that can be used to store the key paramaters about the tested run 
"""
#==============================================================================

__title__ = "Mapping Class"
__author__ = "Arden Burrell"
__version__ = "1.0 (13.02.2018)"
__email__ = "arden.burrell@gmail.com"

#==============================================================================
class mapclass(object):
	"""Class of object that contains infomation about the run
	Takes:
		The opened runs_info csv file 
		The run number
	Returns:
		
	"""
	def __init__(self, region, plotpath, pshow=True):
		# ========== the variables defined by the rundet ========== 

		# Useful infomation
		self.region = region
		aus_names = ["AUS", "Australia", "aus"]
		# ===== Using the run infomation, determine the bounds for a map =====
		if region == ["GLOBAL", "global"]:
			self.bounds = [-180, 180, 90, -90]
		elif region in aus_names:
			self.bounds = [112.0, 156.25, -44.5, -10]
		elif region == "MONG":
			self.bounds = [85.0, 120.0, 52.0, 40.0]
		else:
			Warning("The region code is unknown, unable to set bounds")
			self.bounds = None
		
		self.pshow  = pshow # show the plot after saving?
		
		# ========== Set the blank variables ========== 
		self.desc     = None # The variable being mapped
		self.var      = None # The variable being mapped
		self.cmap     = None # Colormap set later
		self.cmin     = None # the min of the colormap
		self.cmax     = None # the max of the colormap
		self.cZero    = None # the zero point of the colormap
		self.column   = None # the column to be mapped
		self.norm     = None
		self.ticks    = None # The ticks on the colorbar
		self.ticknm   = None # The tick names on the colorbar
		self.cbounds  = None
		self.plotpath = plotpath
		# ========== Set the Histogram lables ==========
		self.cblabel  = None # add label to colorbar
		self.cbsci    = None # The scientific notation power range
		self.dpi      = None # The DPI of the output figures
		self.national = True # Inculde national boundaries
		self.borders  = False # Figure borders
		self.set_x    = 3.1 # The spacing of the colorbar lables
		self.xlocator = None # mticker.FixedLocator(np.arange(80.0, 125.0, 5.0))
		self.ylocator = None # mticker.FixedLocator(np.arange(56.0, 40.0, -2.0))
		self.gridalp  = 0.5  # The alpha of the grid
		self.font     = None # Font family name
		self.text     = None # the text(usually an "a)} or equ )  

		
		
		# ========== Set the plot defualts ==========
		self.maskcol  = 'dimgrey'
		self.Oceancol = 'w'
		# Extend the colorbar to cover any data outside the cmap range
		self.extend    = "both" 
		self.spacing   = 'uniform'
		self.fontsize  = 11
		self.latsize   = 16

		


