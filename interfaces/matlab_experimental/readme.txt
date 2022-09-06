This new Matlab Toolbox for Cantera changs the Matlab interface to the modern Matlab structure and syntaxes for OOP. It replaced the MEX interface with direct function calling from Canter Clib.

Installation guide:

1) Install Matlab (any release newer than R2008a). 
2) Compile Cantera from Source and install in your Python environment, as directed in this link. 
https://cantera.org/install/compiling-install.html. Be sure to set `matlab_path` option in Scons.
3) Both the new and legacy versions of the Matlab toolbox will be installed. However, only one can be used at a time. Activating both will result in a conflict in Matlab. 
4) For first time users, launch Matlab, then navigate to [/path/to/cantera/installation/matlab] using "Browse for Folder".
5) In the Command Window, run `ctpath('install')` to add search path to the Cantera Matlab Toolbox. 
6) The user can then set the desired toolbox version by running the 'SetMatlabToolboxVersion' command. Input can be 'legacy' or 'new'.
7) The user must restart Matlab after switching between the new and legacy interfaces. 
8) To use the new Cantera interface, first set toolbox version to 'new' and run `LoadCantera` command. 
9) To stop using the new Cantera interface, run the following commands:
`clear all`
`close all`
`cleanup`
`UnloadCantera`
10) To use the legacy Cantera interfaces, simply set toolbox version to 'legacy'.
11) TO remove the search path to the Cantera Matlab Toolbox, run `ctpath('uninstall')` command in Matlab. This does NOT uninstall the Toolbox, but is a necessary step before uninstalling Cantera using SCons as SCons does not clean up the search paths for Matlab. 
