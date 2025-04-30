%MPM_ASmakekm : wave number matrix km
% 
% km = MPM_ASmakekm(k,coord,M,bz)
% 
% Calculates the modal wave number matrix km.
% 
% Input parameters: 
% k : free space wave number
% coord : (z, radius) of the tube
% M : number of modes
% bz : zeros of Bessel function J1
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
function km = MMM_ASmakekm(k,coord,M,bz)
R = coord(2);
gamma_m = bz(1:M);
gmR = (gamma_m./R);
gmRm = gmR(:,ones(1,length(k)));
km = k(ones(M,1),:);
km = conj(sqrt(km.^2-gmRm.^2));
