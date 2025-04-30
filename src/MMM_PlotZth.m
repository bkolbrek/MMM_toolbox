%MMM_PlotZth: plots normalized throat impedance
%
% MMM_PlotZth(data, figNo)
%
% Displays a plot of the normalized throat impedance.
%
% Input parameters:
%   data : MMM data struct
%   figNo : the figure number (optional)
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

function MMM_PlotZth(data, figNo, normalize)
if nargin < 3
	normalize = true;
end
if nargin < 2
    figure();
else
    figure(figNo);
end

Z00n = data.Z00;
if normalize
	Z00n = data.St/(data.rho*data.c)*Z00n;
end

semilogx(data.fvec, real(Z00n), 'k', data.fvec, imag(Z00n), 'r');
xlim(data.fvec([1,end]));
grid;
title('Horn throat impedance');
ylabel('Normalized acoustic Z');
xlabel('Hz');