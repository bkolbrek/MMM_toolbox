%MMM_ASplotHorn: plots the horn profile
%
% MMM_ASplotHorn(data, figNo)
%
% Plots the initial and discretized horn profiles.
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

function MMM_ASplotHorn(data, figNo)
if nargin < 2
    figure();
else
    figure(figNo);
end
ym = max(data.rawCoords(:,2));
plot(data.steppedCoords(:,1),data.steppedCoords(:,2), 'b', data.steppedCoords(:,1),-data.steppedCoords(:,2), 'b',...
    data.rawCoords(:,1),data.rawCoords(:,2), 'k', data.rawCoords(:,1),-data.rawCoords(:,2), 'k')
ylim([-ym, ym]*1.1)
axis equal
xlabel('z axis [m]');
ylabel('Radius [m]');
title('Horn profile');
drawnow;