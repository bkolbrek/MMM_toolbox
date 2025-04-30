%MMM_makebigfmat : calculate F matrices throughout a horn
% 
% bigF = MMM_makebigfmat(N, coords, modeinfo, ffunc)
% 
% Calculates the scattering matrices F  for all discontinuities  in the 
% horn. bigF(:,:,i) is the matrix F at position i in the horn.
%  
% Input parameters: 
% N : number of modes (in each direction for rectangular horns)
% coords : horn coordinates
% modeinfo : required info to calculate the matrix, for instance
%	eigenvalues.
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
function bigF = MMM_makebigfmat(N, coords, modeinfo, ffunc)

Lc = size(coords,1);
bigF = zeros(N,N,Lc);
for iz = 1:Lc-1% (length(coords)-1):-1:1
    L = coords(iz+1,1) - coords(iz,1);
    if (L==0) %propagate across discontinuety
        F = ffunc(N, coords(iz,:), coords(iz+1,:), modeinfo); 
        bigF(:,:,iz)=F;
    end
end
