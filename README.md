# PO_database
%=============================
% Introduction
%=============================
This database accompanies a study[1] performed by Luke Bury, Jay McMahon, and Martin Lo. Please see the paper for a full explanation of the results and further acknowledgements. The Matlab[2] script "main_readPlotPOs.m" is provided to visualize the database. Please see the section below for instructions.

%=============================
% Contents - bin
%=============================
The "bin" folder contains spare functions

%=============================
% Contents - data
%=============================
The "data" folder contains .txt files composed of initial conditions for families of periodic orbits in the Jupiter-Europa system in the CR3BP. In each .txt file, the data is arranged in a csv format. Each row represents a family member. The first 6 columns are the classical 6 state variable (x,y,z,x-dot,y-dot,z-dot). The 7th column is the time period for the periodic orbit. The 8th column is the Jacobi constant of the periodic orbit. The 9th and 10th columns contain the two stability indices of the orbit. The 11th and 12th column contain the alpha and beta values that are used for locating bifurcations with the Broucke stability diagram[3]. The 13th column is a flag that indicates if the orbit impacts the surface of Europa (a value of 1 represents an impact, while a 0 represents no impact). The 14th column is the propagation error - see[1] for an explanation of propagation error.

%=============================
% Contents - main_readPlotPOs.m
%=============================
main_readPlotPOs.m is provided as a simple, high-level script for visualizing the data. Use of the script simply setting a few plotting preferences as "true" or "false," selecting the desired family to be plotted, and selecting certain members from that family to plot.

%=============================
% Instructions for using main_readPlotPOs.m
%=============================
1) In section (2), there are two parameters you can turn on or off with "true" or "false."
	-plot_JCVsTp: Setting this parameter to "true" will produce a simple plot where each member of the family is represented on a plot of Jacobi constant vs time period.
	
	-plot_family: Setting this parameter to "true" will produce a 3D plot of certain members of the family. These members will be propagated for one time period and plotted.

2) In section (3-1), choose a font size for plot labels and also an integration tolerance for use with ode

3) In section (3-2), un-comment the family you wish you plot results from

4) In section (4-3-1), set the variable "plot_PO_indices" equal to a vector of indicies of POs you wish to plot from the database. One option is to manually select the indices you wish to plot. This can be done in (4-3-1-1). The alternative is to plot 'n' members equally spaced along the JC-Tp curve of the family. This 'n' value can be set as the variable "n_plotPOs" in (4-3-1-2). 

If you are using the method of (4-3-1-2) to choose plotting indices, be sure to comment-out (4-3-1-1), and vice versa. 

5) Run the script.

%=============================
% References
%=============================
[1] L. Bury, J. McMahon, and M. Lo, “A Study of Periodic Orbits near Europa,” In Progress, 2021.
[2] The Math Works Inc., “MATLAB R2021a”
[3] R. Broucke. Stability of Periodic Orbits in the Elliptic, Restricted Three-Body Prob- lem. AIAA Journal, 7(6):1003–1009, 1969.



