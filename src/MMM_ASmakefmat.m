%MPM_makefmat : Make scattering matrix F
% 
% F = MPM_ASmakefmat(N,coord1,coord2,bz);
% 
% Calculates the coupling matrix F, used to propagate modes across
% a discontinuity (axisymmetric case).
% 
% Returns the inverse of the matrix V if R1>R2
% 
% Input parameters: 
% N : number of modes
% coord1 : (z, radius) of tube 1
% coord2 : (z, radius) of tube 2
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

function F = MMM_ASmakefmat(N,coord1,coord2,bz)
R1 = coord1(2);
R2 = coord2(2);
beta = R1/R2;

if (beta > 1)
    beta = 1/beta;
elseif (beta==1)
    F = eye(N);
    return;
end

F = zeros(N);
    
gamma_m =bz(:,ones(1,N)).';
gamma_n =bz(:,ones(1,N));
Fm = 2*beta*gamma_m.*besselj(1,beta*gamma_m)./ besselj(0,gamma_m);
F = Fm./(beta^2*gamma_m.^2 - gamma_n.^2);

F(1,1) = 1;

