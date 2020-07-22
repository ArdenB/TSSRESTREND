"""
CHeck an see if a folder exists and if not makes directory
"""
#==============================================================================

__title__ = "python mkdir"
__author__ = "Arden Burrell"
__version__ = "1.0 (17.02.2018)"
__email__ = "arden.burrell@gmail.com"

#==============================================================================
import os
#==============================================================================

def pymkdir(path):
	"""Makes a folder if it doens't already exist
	args:
		path
	"""
	if not os.path.exists(path):
		print(path)
		os.makedirs(path)

#==============================================================================