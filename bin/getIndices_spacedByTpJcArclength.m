function [targetIndices] = getIndices_spacedByTpJcArclength(Tps, JCs, n_desired)
%%% Description
% Databases of periodic orbit families are hard to space-out correctly.
% Families can be non-linear and even non-unique in both the elements of Tp
% and JC. The purpose of this code is to scan through a database
% (represented by vectors of the time periods and jacobi constants) and
% determine indices of 'N' number of periodic orbits which are near-evenly
% space along the curve of Tp vs JC. 
%
% This is useful for trying to plot families from databases, and helps to
% circumvent the problem that arrises naturally from populating some
% regions of the family much more densely as continuation algorithms are
% used through sensitive regions that require small steps to advance.
% 
% ------------------------------------------------------------------------
%%% Inputs
% Tps       - [nx1 or 1xn] vector of time periods
% JCs       - [nx1 or 1xn] vector of Jacobi constants
% n_desired - [scalar] number of desired indices in output 
% ------------------------------------------------------------------------
%%% Outputs
% targetIndices - [nx1] vector of indices for input vectors. So for
%                  example, to plot evenly space members from the database,
%                  you can say: PO_data(targetIndices,c_x:c_z)
% ------------------------------------------------------------------------
% Created: 2/11/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Check that dimensions & lengths of two input vectors are the same
% -------------------------------------------------
if ~isequal(size(Tps), size(JCs))
    warning('Sizes or dimensions of inputs don''t match')
    return
end

% -------------------------------------------------
%%% Normalize the span of each input so that neither outweigh the other
% -------------------------------------------------
Tps = Tps - min(Tps);
Tps = Tps ./ max(Tps);

JCs = JCs - min(JCs);
JCs = JCs ./ max(JCs);

% -------------------------------------------------
%%% Calculate the arclength of each new point on the [Tp, Jc] curve
% -------------------------------------------------
%%% Preallocate
numericalArclengths           = NaN(length(Tps)-1,1);
progressingNumericalArclength = NaN(length(Tps)-1,1);

for kk = 2:length(Tps)
    %%% Length of each segment
    numericalArclengths(kk-1)           = ((Tps(kk) - Tps(kk-1))^2 + (JCs(kk) - JCs(kk-1))^2)^(1/2);
    
    %%% New total length at the current segment
    progressingNumericalArclength(kk-1) = sum(numericalArclengths(1:kk-1));
end

% -------------------------------------------------
%%% Calculate specific arclengths that are nearly evenly spaced along the
%%% [Tp, Jc] curve based on its total length
% -------------------------------------------------
arclengthTargets = linspace(0, progressingNumericalArclength(end), n_desired);

% -------------------------------------------------
%%% For each target arclength, find the index that most nearly matches it
% -------------------------------------------------
targetIndices = NaN(n_desired, 1);
for kk = 1:n_desired
    searchVector = abs(progressingNumericalArclength - arclengthTargets(kk));
    temp = find(searchVector == min(searchVector));
    targetIndices(kk) = temp(1);
end

%%% Re-adjust output vector of indices so that it starts on the first index
%%% of the input vectors and ends on the last index
targetIndices = targetIndices + 1;
targetIndices(1) = 1;

if length(targetIndices) ~= length(unique(targetIndices))
	warning('Database not sufficiently dense for desired number of outputs')
end

end % function