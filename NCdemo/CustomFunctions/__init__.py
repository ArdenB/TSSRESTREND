import sys

if sys.version.startswith("2.7"):
	from MakeDirectory import pymkdir
	from WriteMetadata import writemetadata
	from GitMetadata   import gitmetadata
	from MapMaker      import mapmaker
	from MapClass      import mapclass
	from replacehex    import ReplaceHexColor
elif sys.version.startswith("3."):
	from CustomFunctions.MakeDirectory import pymkdir
	from CustomFunctions.WriteMetadata import writemetadata
	from CustomFunctions.GitMetadata   import gitmetadata
	from CustomFunctions.MapMaker      import mapmaker
	from CustomFunctions.MapClass      import mapclass
	from CustomFunctions.replacehex    import ReplaceHexColor
else:
	print("Unknow Version of python")
	raise ImportError

__all__ = [ "pymkdir", "writemetadata"]