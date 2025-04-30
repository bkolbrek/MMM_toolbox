%MMM_ASpolarMap: plots polar map and directivity index
%
% MMM_ASpolarMap(data, angles, normalize, plotDI, figNo)
%
% Displays the polar map and directivity index (optional). The polar map
% can be normalized. 
%
% Input parameters:
%   data : MMM data struct
%   angles : angles of the field points
%   normalize : normalize the polar map
%   plotDI : include a plot of the directivity index
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

function MMM_ASpolarMap(data, angles, normalize, plotDI, figNo)
if nargin < 5
    figure();
else
    figure(figNo);
end
clf;
if nargin < 4
    plotDI = true;
end
if nargin < 3
    normalize = true;
end

Lp = 20*log10(abs(data.pRad));
if normalize
    Lp = Lp - Lp(1,:);
end

if plotDI
    subplot(2,1,1);
end
contourf(data.fvec, angles, Lp, 15);
set(gca, 'xscale', 'log');
xlim(data.fvec([1,end]));
ylabel('Degrees');
xlabel('Hz');
t = 'Polar map';
if normalize
    t = [t ' (normalized)'];
end
title(t);
colorbar;

if plotDI    
    % calculate DI
    data = MMM_ASgetDI(data, angles);
    figure(figNo);
    subplot(2,1,2);
    semilogx(data.fvec, data.DI);
    xlim(data.fvec([1,end]));
    xlabel('Hz');
    ylabel('dB');
    title('Directivity index');
    grid;
end
