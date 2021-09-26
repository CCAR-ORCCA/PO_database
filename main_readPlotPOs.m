% ========================================================================
%%% Description
% ========================================================================
% This script is used to read families of periodic orbits from a database
% and plot desired information.

% Created: July 19, 2021
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% (1) Initialization
% ========================================================================
clear
clc
close all
addpath('bin')
ticWhole = tic;

% ========================================================================
%%% (2) Run-Switches
% ========================================================================
plot_JCVsTp = true; % Plot Jacobi constant vs Time period
plot_family = true; % Plot 3D visualizations of family members 

% ========================================================================
%%% (3) Setup
% ========================================================================
% -------------------------------------------------
%%% (3-1) Quick Options
% -------------------------------------------------
%%% Font size for plot axis labels
fontSize = 26;

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% (3-2) Choose PO family
% -------------------------------------------------
% family = 'Jupiter_Europa.CR3BP.L2_Lyapunov.txt';
% family = 'Jupiter_Europa.CR3BP.L2_L_1T.txt'; % south halo family
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_3P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2.txt'; % Butterfly
% % % % % % family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P4.txt';
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P3.txt'; % the end of this is the 2P2 bifurcation from LoPO_2P2_1T (so this is also LoPO_2P2_1T_2P2)
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P2.txt';
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_2T.txt'; % Full Axial family (Also L2_V_T)
% family = 'Jupiter_Europa.CR3BP.L2_L_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_Vertical.txt';

% family = 'Jupiter_Europa.CR3BP.DRO.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P4_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P4_1P4_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P3_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P3_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P3_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P4_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P4_2T.txt'; 

% family = 'Jupiter_Europa.CR3BP.DPO.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P4_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_2P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_3P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_3P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_3P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_4P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_4P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_4P2.txt';  % THIS IS DPO_3P2 FROM THE OTHER SIDE
% family = 'Jupiter_Europa.CR3BP.DPO_5P2.txt';  % THIS IS LoPO_3P2 FROM THE OTHER SIDE
% family = 'Jupiter_Europa.CR3BP.DPO_2T.txt'; 

% family = 'Jupiter_Europa.CR3BP.LoPO.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_2P2.txt'; 
family = 'Jupiter_Europa.CR3BP.LoPO_2P2_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_3P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_3P4.txt'; 

% family = 'Jupiter_Europa.CR3BP.Hg1.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_2T.txt'; 

% family = 'Jupiter_Europa.CR3BP.Hg2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg2_1P2.txt';
% family = 'Jupiter_Europa.CR3BP.Hg2_2P2.txt';
% family = 'Jupiter_Europa.CR3BP.Hg2_1P3.txt';
% family = 'Jupiter_Europa.CR3BP.Hg2_2T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg2_3T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg2_5T.txt';

% family = 'Jupiter_Europa.CR3BP.Se7.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_2T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_3T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_5P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_6P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_5T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_6T.txt'; 

% -------------------------------------------------
%%% (3-3) Get system parameters
% -------------------------------------------------
%%% Get primary and secondary as structs, and get normalizing parameters
[primary, secondary, rNorm, tNorm, vNorm] = get_systemParameters(family);

%%% prms for integration
prms.u  = secondary.MR;

if contains(family, '.CR3BP.')
    prms.n = 1;
else
    warning('A perturbed CR3BP model will require a new mean motion value')
    return
end
% -------------------------------------------------
%%% (3-4) Load the family of POs
% -------------------------------------------------
PO_datafile = [pwd,'/data/' family];

%%% Load the data file
PO_data = dlmread(PO_datafile,',',1,0);

%%% Grab header line
fid = fopen(PO_datafile, 'rt');  %the 't' is important!
header = fgetl(fid);
fclose(fid);

%%% Number of ICs
n_POs = size(PO_data,1);

%%% Indices for plotting purposes
PO_indices = linspace(1, n_POs, n_POs);

% --------------------------
% (3-4-1) Create column specifiers
% --------------------------
PO_header_2020 = 'x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag,error';

if contains(header, PO_header_2020)
    c_x0 = 1;   c_y0 = 2;   c_z0 = 3;
    c_xd0 = 4;  x_yd0 = 5;  c_zd0 = 6;
    c_Tp = 7;   c_JC = 8;   c_S1 = 9;   c_S2 = 10;
    c_alpha = 11;   c_beta = 12;    c_impactFlag = 13;
    c_error = 14;
end

% ========================================================================
%%% (4) Plot Data
% ========================================================================
% -------------------------------------------------
%%% (4-1) Setting color data
% -------------------------------------------------
colors.black      = [0 0 0]./255;
colors.blue       = [50  100 200]./255;
colors.blue2      = [47 130 255]./255;
colors.brown      = [140 86 75]./255;
colors.cyan       = [0 255 255]./255;
colors.drkgrey    = [57 59 59]./255;
colors.drkgrn     = [15 87 20]./255;
colors.drkgrn2    = [0 158 26]./255;
colors.drkred     = [150 40 32]./255;
colors.grey       = [127 127 127]./255;
colors.grn        = [44 160 44]./255;
colors.grn2       = [44 235 44]./255;
colors.khaki      = [219 219 141]./255;
colors.ltblue     = [125 216 255]./255;
colors.ltbrown    = [196 156 148]./255;
colors.ltgrey     = [199 199 199]./255;
colors.ltgrn      = [152 223 138]./255;
colors.ltmag      = [229 194 237]./255;
colors.ltorange   = [255 187 120]./255;
colors.ltpink     = [247 182 210]./255;
colors.ltpurp     = [197 176 213]./255;
colors.ltred      = [255 152 150]./255;
colors.ltturq     = [158 218 229]./255;
colors.mag        = [255  0 255]./255;
colors.orange     = [255 127 14]./255;
colors.pink       = [227 119 194]./255;
colors.pink2      = [255 7 178]./255;
colors.purp       = [148 103 189]./255;
colors.purp2      = [146, 0, 250]./255;
colors.red        = [209   0   0]./255;
colors.red2       = [255 91 69]./255;
colors.red3       = [255 53 17]./255;
colors.turq       = [23 190 207]./255;
colors.ylwgrn     = [188 189 34]./255;
colors.ylw        = [247 202   0]./255;
colors.white      = [255 255 255]./255;

% -------------------------------------------------
%%% (4-2) Create plot of Jacobi constant vs time period
% -------------------------------------------------
if plot_JCVsTp
    %%% Plot Jacobi constant and Tp
    figure; hold all
    plot3(PO_data(:, c_Tp), PO_data(:, c_JC), PO_indices,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue)
    
    %%% Make background of plot white
    set(gcf,'color','white')

    %%% Turn on the grid
    grid on

    %%% Set labels, font size, and interpreter
    xlabel('$T_P$','FontName','Times New Roman','Fontsize',fontSize,'Interpreter','LaTex')
    ylabel('Jacobi Constant','FontName','Times New Roman','Fontsize',fontSize,'Interpreter','LaTex')
end


% -------------------------------------------------
%%% (4-3) Create 3D plot of periodic orbits
% -------------------------------------------------
if plot_family
    % --------------------------
    % (4-3-1) Choose indices to plot from database
    % --------------------------
    
    %%% (4-3-1-1) To set specific indices from the database to plot:
%     plot_PO_indices = 1:n_POs;
%     plot_PO_indices = 1;
%     plot_PO_indices = n_POs;
%     plot_PO_indices = [1, n_POs];
    
    %%% (4-3-1-2) To plot 'n' POs which are evenly spaced along the Tp-JC curve of the family
    n_plotPOs = 10;
    
    plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), n_plotPOs);
    
    
    % --------------------------
    % (4-3-2) Propagate the POs corresponding to the selected indices and plot
    % --------------------------
    %%% Initialize STM as column vector to be added with state
    stm0_colVec = reshape(eye(6),36,1);
    
    %%% Loop, integrate, and plot
    figure; hold all
    axis equal
    for index = plot_PO_indices'
        [T_PO, X_PO] = ode113(@integrator_CR3BP_STM, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options, prms);
        
        plot3(X_PO(:,1), X_PO(:,2), X_PO(:,3), 'linewidth', 1, 'color', colors.blue2) % 0.1961    0.3922    0.7843
    end
    
    %%% Set labels, font size, and interpreter  
    xlabel('$x_n$','FontName','Times New Roman','Fontsize',fontSize,'Interpreter','LaTex')
    ylabel('$y_n$','FontName','Times New Roman','Fontsize',fontSize,'Interpreter','LaTex')
    zlabel('$z_n$','FontName','Times New Roman','Fontsize',fontSize,'Interpreter','LaTex')

    %%% Make background of plot white
    set(gcf,'color','white')

    %%% Turn on the grid
    grid on
end % plot_family


% ========================================================================
%%% (5) Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\nElapsed time: %1.4f seconds\n',tocWhole)











