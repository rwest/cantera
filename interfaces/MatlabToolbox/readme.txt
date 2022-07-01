This new Matlab Toolbox for Cantera changs the Matlab interface to the modern Matlab structure and syntaxes for OOP. It replaced the MEX interface with direct function calling from Canter Clib.

Installation guide:

1) Install Matlab (any version newer than R2008a). The toolbox have been tested for R2020a, R2020b, R2021a, and R2021b. 
2) Compile Cantera from Source and install in your Python environment, as directed in this link. 
https://cantera.org/install/compiling-install.html
3) Both the new and legacy versions of the Matlab toolbox will be installed. However, only one can be used at a time. Activating both will result in a conflict in Matlab. 
4) For first time users, launch Matlab, then navigate to [/path/to/cantera/installation/matlab/MatlabToolbox] using "Browse for Folder".
5) The user can set the desired toolbox version by running either `ctpath_new.m` (for the new version) or `ctpath_legacy.m` (for the old version), then save the search path as the default search path in Matlab using "Set Path". This process is reversible, and the user can run this step to switch between the two versions at any time.  
6) To use the new Cantera interface, launch Matlab and run `LoadCantera` command. 
7) To stop using the new Cantera interface, run the following commands:
`clear all`
`close all`
`cleanup`
`UnloadCantera`
