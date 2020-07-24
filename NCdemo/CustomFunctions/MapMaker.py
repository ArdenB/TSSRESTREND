# -*- coding: utf-8 -*-
"""
Function to make maps


"""
#==============================================================================

__title__ = "Map Maker"
__author__ = "Arden Burrell"
__version__ = "1.2 (11.03.2018)"
__email__ = "arden.burrell@gmail.com"

#==============================================================================

# Import packages
import numpy as np
import pandas as pd
import sys	
import ipdb
import xarray as xr
# import datetime as dt

# Mapping packages
import cartopy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cpf
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl 
import matplotlib.font_manager as font_manager
# Make sure it waits to render plots
mpl.interactive(False)

import matplotlib.pyplot as plt
import warnings as warn

# if sys.version.startswith("2.7"):
# 	from  .. import CoreFunctions as cf 	
# elif sys.version.startswith("3."):
# 	from  .. import CoreFunctions as cf 	

#==============================================================================

#==============================================================================
def mapmaker(da, mapdet):
	"""Function to Build some maps"""

	plt.rcParams.update({'figure.subplot.right' : 0.82 })
	plt.rcParams.update({'figure.subplot.left' : 0.02 })
	plt.rcParams.update({'figure.subplot.top' : 0.95 })
	plt.rcParams.update({'figure.subplot.bottom' : 0.05 })
	# plt.rcParams.update({'figure.subplot.right' : 0.80 })
	# plt.rcParams.update({'figure.subplot.left' : 0.10 })
	
	# # ========== Deal with the fonts ==========
	# if not mapdet.font is None:
	# 	font_dirs = ['/mnt/c/Windows/Fonts', ]
	# 	font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
	# 	font_list = font_manager.createFontList(font_files)
	# 	# breakpoint()
	# 	# font_list = font_manager.FontManager.addfont(font_files)
	# 	font_manager.fontManager.ttflist.extend(font_list)
	# 	# breakpoint()
	# 	font = {'family' : mapdet.font,
	# 			'weight' : 'bold', #,
	# 			'size'   : mapdet.latsize}
	# else:
	font = {'weight' : 'bold', #,
			'size'   : mapdet.latsize}

	matplotlib.rc('font', **font)
	plt.rcParams.update({'axes.titleweight':"bold", 'axes.titlesize':mapdet.latsize})

	aus_names = ["AUS", "Australia"]



	if mapdet.region in aus_names:
		fw =12.8
		fh =9.6
		# fig, ax = plt.subplots(1, 1, figsize=(12,9),
		# 	subplot_kw={'projection': ccrs.PlateCarree()}, 
		# 	num=("Map of %s" % mapdet.var), dpi=mapdet.dpi)
	if mapdet.region in "MONG":
		fw =18
		fh=6.5
		# fig, ax = plt.subplots(1, 1, figsize=(18,6.5),
		# 	subplot_kw={'projection': ccrs.PlateCarree()}, 
		# 	num=("Map of %s" % mapdet.var), dpi=mapdet.dpi)

	else:
		fw =18
		fh = 8
		# fig, ax = plt.subplots(1, 1, figsize=(18,8),
		# 	subplot_kw={'projection': ccrs.PlateCarree()}, 
		# 	num=("Map of %s" % mapdet.var), dpi=mapdet.dpi)

	fig= plt.figure(
		figsize=(fw, fh),
		num=("Map of %s" % mapdet.var), dpi=mapdet.dpi)
	ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

	if mapdet.region == "GLOBAL":
		ax.set_global()
	else:
		ax.set_extent(mapdet.bounds, crs=ccrs.PlateCarree(),)
	
	# ipdb.set_trace()
	# ========== plot the image ==========
	# if mapdet.region == "MONG":
	# 	origin="upper"
	# else:
	im = da.isel(time=0).plot.imshow(
		cmap=mapdet.cmap, 
		norm=mapdet.norm, 
		ax  = ax, 
		add_colorbar = False, 
		extend ="neither",
		vmin = mapdet.cmin,
		vmax = mapdet.cmax
		) # added after australia looked lame
	ax.set_title(mapdet.var)
	# re-calculated at each figure resize. 
	posn = ax.get_position()
	cbar_ax = fig.add_axes([posn.x0 + posn.width + 0.005, posn.y0, 0.025, posn.height])


	# ========== Add features to the map ==========
	ax.add_feature(cpf.OCEAN, facecolor="w", alpha=1, zorder=100)
	ax.add_feature(cpf.COASTLINE, zorder=101)
	if mapdet.national:
		ax.add_feature(cpf.BORDERS, linestyle='--', zorder=102)
	ax.add_feature(cpf.LAKES, alpha=0.5, zorder=103)
	ax.add_feature(cpf.RIVERS, zorder=104)
	ax.outline_patch.set_visible(False)
	# ax.gridlines()

	# =========== Set up the axis ==========
	gl = ax.gridlines(
		crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=mapdet.gridalp, 
		linestyle='--', zorder=105)

	gl.xlabels_top = False
	gl.ylabels_right = False

	# gl.xlines = False
	if mapdet.region in aus_names:
		gl.xlocator = mticker.FixedLocator(range(110, 170, 10))
		gl.ylocator = mticker.FixedLocator(range(-10, -60, -10))
	elif mapdet.region == "MONG":
		gl.xlocator = mticker.FixedLocator(np.arange(80.0, 125.0, 5.0))
		gl.ylocator = mticker.FixedLocator(np.arange(56.0, 40.0, -2.0))
	# elif mapdet.region == "Cropped":
	# 	gl.xlocator = mticker.FixedLocator(range(-140, 161, 20))
		# gl.ylocator = mticker.FixedLocator(range(mapdet.crop[0], mapdet.crop[1], -10))
	if not mapdet.xlocator is None:
		gl.xlocator = mticker.FixedLocator(mapdet.xlocator)
	else:pass

	if not mapdet.ylocator is None:
		gl.ylocator = mticker.FixedLocator(mapdet.ylocator)
	else:pass

	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER
	# gl.xlabel_style = {'size': 15, 'color': 'gray'}
	# gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
	
	# ========== Add an autoresizing colorbar ========== 
	def resize_colobar(event):
		plt.draw()

		posn = ax.get_position()
		cbar_ax.set_position([posn.x0 + posn.width + 0.01, posn.y0, 0.024, posn.height])

	fig.canvas.mpl_connect('resize_event', resize_colobar)

	posn = ax.get_position()
	cbar_ax.set_position([posn.x0 + posn.width + 0.010, posn.y0*1.1, 0.024, posn.height*0.97]) 

	# ========== save the plot ==========
	if not mapdet.text is None:
		ax.text(-0.04, 0.96, mapdet.text, transform=ax.transAxes, 
					size=mapdet.latsize, weight='bold')#, zorder=106)
	cb = plt.colorbar(
		im, 
		cax        = cbar_ax, 
		extend     = mapdet.extend, 
		norm       = mapdet.norm,
		ticks      = mapdet.ticks, 
		spacing    = mapdet.spacing,
		boundaries = mapdet.cbounds
		)

	if not (mapdet.ticknm is None):
		cb.ax.set_yticklabels(mapdet.ticknm) 
		# cb.ax.xaxis.set_tick_params(pad=25)

		for t in cb.ax.get_yticklabels():
			t.set_fontsize(mapdet.fontsize)
		# ipdb.set_trace()
	else:
		# Change the horrixontal allignment of cb ticks to right
		# cb.ax.xaxis.set_tick_params(pad=25)
		for t in cb.ax.get_yticklabels():
			t.set_horizontalalignment('right')   
			# sx = 3.1
			t.set_x(mapdet.set_x)
			# t.set_x(3.0)
			t.set_fontsize(13)
	
	if not (mapdet.cblabel is None):
		cbar_ax.set_ylabel(mapdet.cblabel, rotation=90, weight='bold')

	if mapdet.borders:
		#Adding the borders
		for bord in ['top', 'right', 'bottom', 'left']:

			ax.spines[bord].set_visible(True)
			ax.spines[bord].set_zorder(20)
			ax.spines[bord].set_color('k')
			ax.spines[bord].set_linewidth(2.0)

	if mapdet.dpi is None:
		dpi = fig.dpi
	else:
		dpi = mapdet.dpi
	plt.draw()

	# ========== save the plot ==========
	if not (mapdet.plotpath is None):
		fname = "%sMapof%s" % (
			mapdet.plotpath, mapdet.desc)
		for ext in [".png", ".pdf"]:
			# Make a pdf version
			plt.savefig(fname+ext, dpi=fig.dpi)
			
		plotinfo = "PLOT INFO: Plot of %s made using %s:v.%s" % (mapdet.var, __title__, __version__)

	plt.show()
	# ipdb.set_trace()
	# Reset the the plt paramters to the defualts
	plt.rcParams.update(plt.rcParamsDefault)
	# return the infomation
	return plotinfo, fname

