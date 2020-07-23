# -*- coding: utf-8 -*-

""" 
Create my own style of colorbar

"""

__title__ = "Replace Hex Colors"
__author__ = "Arden Burrell"
__version__ = "1.1(19.02.2018)"
__email__ = "arden.burrell@gmail.com"

#==============================================================================
# modules
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mpc
# import webcolors as wc
import ipdb

#==============================================================================

def ReplaceHexColor(cmaplist, cmin, cmax, color = '#FFFFFF', ticks=None, zeroR = None):
	""" 
	Take a list of hex values, insert white in as a color, define some bounds
	args:
		cmapList

	"""
	# ========== Setup the parpmaters ==========
	# Get the size
	cmapsize  = int(len(cmaplist))
	# get the number of ticks
	tickinter = (cmax-cmin)/(float(cmapsize))

	# ========== add the color to the list ==========
	cmaplist.insert(int(cmapsize/2), color)

	# =========== Remove the weak colors ========== 
	# del cmaplist[cmapsize/2+1]
	# del cmaplist[cmapsize/2-1]
	# cmaplist.insert(0, cmaplist[0])
	# cmaplist.insert(-1, cmaplist[-1])
	# pdb.set_trace()
	if zeroR is None:
		if cmax < 0.01:
			zerob = 0.000001
		else:
			zerob = 0.0001
	else:
		zerob = zeroR
	
	# Set the tick values
	try:
		if ticks is None:
			ticks 	  = np.round(np.arange(cmin, (cmax+0.00001), tickinter), decimals=6)
			spacing = 'proportional'
			bounds = np.hstack(
				(
				ticks[:int(cmapsize/2)], 
				np.array([-zerob, zerob]), 
				ticks[int(cmapsize/2+1):]
				)
			)
		else:
			spacing  = 'uniform'
			bounds = np.hstack(
				(
				ticks[:int(cmapsize/2)], 
				np.array([-zerob, zerob]), 
				ticks[int(cmapsize/2+1):]
				)
			)
			ticks = np.hstack(
				(ticks[:int(cmapsize/2)], 
				np.array([-zerob,ticks[int(cmapsize/2)], zerob]), 
				ticks[int(cmapsize/2+1):]))
	
	except Exception as e:
		print(str(e))
		ipdb.set_trace()
		

	# ========== create the colormap ==========
	cmap   = mpl.colors.ListedColormap(cmaplist)
	cmap.set_over(cmaplist[-1])
	cmap.set_under(cmaplist[0])
	
	# ========== Define the bounds and norm ==========
	
	
	norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)
	
	
	# ========== Return cmap and infomation ==========
	return cmap, norm, ticks, bounds, spacing