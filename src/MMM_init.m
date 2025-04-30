%MMM_init: initialize MMM data struct
%
% data = MMM_init(fvec, rho, c, nModes, eigenValues, coords)
%
% Initializes the data struct that contains the calculation data.
% If the input data contains coordinates, stepped coordinates are generated
% and the F-matrix calculated.
%
% Input parameters:
%   fvec: frequency vector [Hz]
%   nModes : maximum number of modes
%   coords : horn coordinates (not stepped)
%	geometry : type of geometry:
%		'axi' : axisymmetric
%		'rect' : rectangular
%   rho : air density (optional, default = 1.205)
%   c : sound speed (optional, default = 344)
%
%  --------------------------------------------------- -------------------------------------------
%    This file is part of the Mode Matching Method (MMM) Toolbox by Bjørn Kolbrek.
%    Copyright (C) 2012-2025 by Bjørn Kolbrek
%       https://kolbrek.hornspeakersystems.info/
%		https://github.com/bkolbrek/MMM_toolbox
%
%    The MPM Toolbox is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by the Free Software
%    Foundation, either version 2 of the License, or (at your option) any later version.
%
%    The MPM Toolbox is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License along with the
%    MPM Toolbox. If not, see <http://www.gnu.org/licenses/>.
%  --------------------------------------------------- -------------------------------------------
%
function data = MMM_init(fvec, nModes, coords, geometry, rho, c)
if nargin < 5
	rho = 1.205;
end
if nargin < 6
	c = 344;
end
data = struct();

if isempty(coords)
	error('Error: no horn coordinates.')
elseif any(diff(coords(:,1))<0)
	error('Error: Z-coordinate is not monotonically increasing.');
end


data.geometry = lower(geometry);

data.rho = rho;
data.c = c;
data.fvec = fvec;
data.nfreq = length(fvec);
data.k = fvec*2*pi/data.c;
data.nModes = nModes;
data.keepZmatrix = true;
data.nIntegrationPoints = 20; % number of radial integration points for rayleigh integral
data.rawCoords = coords;
data.steppedCoords = MMM_makesteps(coords);
data.modeIndex = [];

if contains(data.geometry, 'axi')
	load('MMM_besselzeros.mat','bz');
	data.eigenValues = bz(1:nModes);
	data.S = pi*data.steppedCoords(:,2).^2;
	data.makekm = @MMM_ASmakekm;
	data.makefmat = @MMM_ASmakefmat;
	data.modeInfo = data.eigenValues;
elseif contains(data.geometry, 'rect')
	data.S = data.steppedCoords(:,2) .* data.steppedCoords(:,3) * 4;
	data.eigenValues = [];
	data = MMM_REaddModeIndex(data);
	data.makekm = @MMM_REmakekm;
	data.makefmat = @MMM_REmakefmat;
	data.modeInfo = data.modeIndex;
end
data.Sm = data.S(end);
data.St = data.S(1); 
data.bigF = MMM_makebigfmat(nModes, data.steppedCoords, data.modeInfo, data.makefmat);
data.Zrad = [];