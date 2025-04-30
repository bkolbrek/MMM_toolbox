%MMM_makesteps(coords)
%
% newcoords = MMM_makesteps(coords)
%
% Given a list of coordinates, creates a stepped approximation of the
% profile. The "step point" is midway between each coordinate point.
%
% NOTE: THIS FUNCTION MUST BE RUN BEFORE ANY HORN SIMULATION,
%       to make sure the coordinates are correctly formatted.
%
% Input parameters:
% coords: list of coordinates. First coordinate must be the z coordinate.
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
function newcoords = MMM_makesteps(coords)
s = size(coords);
N = s(1);
D = s(2);
if D == 3
    s(1) = s(1)*4;
else
    s(1) = s(1)*2;
end
newcoords = zeros(s);
i = 1;
newcoords(i,:) = coords(1,:);
for ih = 1:N-1
    i = i+2;
    newcoords(ih*2,1) = coords(ih,1) + (coords(ih+1,1) - coords(ih,1))/2;
    newcoords(ih*2,2:end) = coords(ih,2:end);
    newcoords(ih*2+1,1) = newcoords(ih*2,1);
    newcoords(ih*2+1,2:end) = coords(ih+1,2:end);
end
i = i+1;
newcoords(i,:) = coords(end,:);
if D == 3
    newcoords = newcoords(1:i,:);
end

