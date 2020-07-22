# -*- coding: utf-8 -*-

""" 
Write Metadata about a saved file or plot

"""

__title__ = "Write Metadata about a save file"
__author__ = "Arden Burrell"
__version__ = "1.0 (18.02.2018)"
__email__ = "arden.burrell@gmail.com"

def writemetadata(fname, infomation):
	"""
	Aggs:
		fname filename that the infomation is being written for 
		infomation: list of infomation 
	"""
	if fname[-1] == ".":
		f = open('%stxt' % fname,'w')
	else:
		f = open('%s.txt' % fname,'w')
	for info in infomation:
		f.write("%s\n" %info)
	f.close()
	