%MMM_ASgeteigenfunctions
%
% phi = MMM_ASgeteigenfunctions(R, rcoords,  eigenValues, normalize)
% 
% Calculates the eigenfunctions for a duct of radius R, at the radiuses
% given in rcoords. 
%
% Input parameters: 
% R : max radius of the duct
% rcoords : radial coordinates where the eigenfunction values are
%      calculated. 
% eigenValues : zeros of Bessel function J1
% normalize : (optional, boolean) normalizes the eigenfunctions to the
%      value at R. Default true.
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
function phi = MMM_ASgeteigenfunctions(R, rcoords,  eigenValues, normalize)
if nargin<4 
    normalize = true;
end
alpha = eigenValues/R;
Nr = length(rcoords);
gamma_n = eigenValues(:,ones(1,Nr))';
phi = (besselj(0, alpha*rcoords))';
if normalize
    norm = besselj(0,gamma_n);
    phi = phi./norm;
end