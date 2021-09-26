function [primary, secondary, rNorm, tNorm, vNorm] = get_systemParameters(familyTag)
%%% Description
%       Assigns structs of body data to a primary and secondary body meant
%       for use in the CR3BP. 
%       
% ------------------------------------------------------------------------
%%% Inputs
% familyTag - [str] A string containing the names of the desired
%                   primary and secondary bodies according to the
%                   convention 'Primary_Secondary'. For example,
%                   familyTag might be
%                   "Jupiter_Europa.CR3BP.L2_Vertical"
% ------------------------------------------------------------------------
%%% Outputs
% primary   - [struct] Struct of data for primary body
% secondary - [struct] Struct of data for secondary body
% rNorm     - [scalar] distance-normalizing constant (km)
% tNorm     - [scalar] time-normalizing constant (sec)
% vNorm     - [scalar] velocity-normalizing constant (km/sec)
% ------------------------------------------------------------------------
% Created: July 19, 2021
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% ------------------------------------------------------------------------
%%% Constants
% ------------------------------------------------------------------------
G     = 6.6726e-20; % km^3 * kg^-1 * s^-2
AU_km = 1.495978707e8; % km

% -------------------------------------------------
%%% Assign Primary
% -------------------------------------------------
if contains(lower(familyTag),'jupiter_')
    %%% Jupiter
    primary.name  = 'jupiter';
    primary.title = 'Jupiter';
    primary.mass  = 1.89819e27; % kg
    primary.u     = primary.mass * G; % km^3 * s^-2
    primary.color = [235,214,173]./255;
%     primary.img   = imread(['bin/textures/jupiterSurfTex.jpg']);
    primary.a     = 778.57e6; % km, https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
    primary.R     = 69911; % km
    primary.J2    = 0.0146956247707; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op / MONTE - Jupiter solution 310
    primary.J4    = -5.913138887463e-04; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op / MONTE - Jupiter solution 310
    primary.J6    = 2.077510523749e-05; % MONTE - Jupiter solution 310

else
    warning('No primary')
    return
end

% -------------------------------------------------
%%% Assign Secondary
% -------------------------------------------------
if contains(lower(familyTag),'_europa')
    secondary.name    = 'europa';
    secondary.title   = 'Europa';
    secondary.color   = [0, 1, 1];
%     secondary.img     = imread(['bin/textures/europaSurfTex.jpg']);
    secondary.a       = 671100; % km
    secondary.R       = 1560.8; % km
    secondary.R_n     = secondary.R / secondary.a;
    secondary.Tp      = 3.551181*86400; % sec
    secondary.meanMot = 2*pi / secondary.Tp; % rad/s
    secondary.MR      = 2.528017528540000e-05; % Value from Mar Vaquero and Javier Roa Vicens at JPL (12/2020)
    secondary.J2      = 4.333968885716003e-04; %(1.938209808166e-04)*sqrt(5); % from MONTE ... sqrt(5) is normalizer (vallado pg 547)
    secondary.mass    = secondary.MR * primary.mass / (1 - secondary.MR); % kg
    secondary.u       = secondary.mass * G; % km^3 / s^2
    secondary.LyapTp  = 3.076475265314829;
else
    warning('No Secondary')
    return
end

% -------------------------------------------------
%%% Assign Normalizing Constants
% -------------------------------------------------
%%% Determine theoretical time period of circular-restricted system
Tp_theoretical = 2*pi*sqrt((secondary.a^3) / (G*(secondary.mass + primary.mass))); % sec

%%% Set distance-normalizing constant
rNorm = secondary.a;         % n <-> km

%%% Set time-normalizing constant
tNorm = Tp_theoretical / (2*pi); % n <-> sec

%%% Set velocity-normalizing constant
vNorm = rNorm / tNorm; % n <-> km/sec


end % function