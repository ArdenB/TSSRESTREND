"""
When a map is made, this script makes a textfile of metadata to accompany it


"""
#==============================================================================

__title__ = "Map Metadata"
__author__ = "Arden Burrell"
__version__ = "1.0 (16.02.2018)"
__email__ = "arden.burrell@gmail.com"

#==============================================================================
# Import packages
import os
import sys
import datetime
import git
import os
import socket
# import ipdb

def gitmetadata():
	"""Create a new history record.

	# To Do, 
		get the __title__, __author__, __version__, __email__
		get the run number
		get the paramter being mapped
		get the map's save name
		get the time stamp
		write our the file
	"""
	# ========== Get the git history record ==========
	time_stamp = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
	exe = sys.executable # get the python version
	args = " ".join(sys.argv) #join all the elements
	repo_dir = os.getcwd()
	try:
		# check the host
		host = socket.gethostbyaddr(socket.gethostname())[0]
		git_hash = git.Repo(repo_dir).heads[0].commit # get the current head
	except git.exc.InvalidGitRepositoryError:
		print('To record the git hash, must run script from top of directory tree in git repo')
		git_hash = 'unknown'
	except TypeError:
		print('To record the git hash, must run script from top of directory tree in git repo')
		git_hash = 'unknown'
		
	entry = """%s: %s %s (Git hash: %s)""" % (time_stamp, exe, args, str(git_hash)[0:7])
	
	return entry