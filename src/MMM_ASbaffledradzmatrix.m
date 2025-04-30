%MMM_ASbaffledradzmatrix : Modal radiation impedance matrix
%
% Zmat = MMM_ASbaffledradzmatrix(k, rho, c, S, maxmodes, bz)
% 
% Calculates the modal radiation impedance matrix for the 
% end of a circular tube terminated in an infinite baffle
% by numerical interation.
% The fundamental (plane wave) mode impedance is calcualted
% by analytical functions.
%
% The matrix Zmat is a square, symmetrical matrix. 
%
% Input parameters:
%   k : wavenumber
%    rho : density of medium
%   c : sound speed in medium
%   S : area of tube opening
%   maxmodes : maximum number of modes calculated
%   bz : zeros of Bessel function J1
%   useHFapprox : (default false): uses an approximation above 2.5x mode
%       cutoff, see PhD thesis for details.
%   progressReport : (default false): prints a progress report.
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

function Zmat = MMM_ASbaffledradzmatrix(k, rho, c, S, maxmodes, bz, useHFapprox, progressReport)
if (nargin < 8)
    progressReport = false;
end
if (nargin < 7)
    useHFapprox = false;
end
Zmat = zeros(maxmodes, maxmodes, length(k));
a = sqrt(S/pi);
kR = k*a;
tol = 10.^(-log(kR)-1);
tol = max(min([0.01 tol]),1e-6);
if length(k) > 1 
    tol = 1e-6;
end

R00 = 1-besselj(1,2*kR)./(kR);
X00 = 2*struveH1(2*kR)./(2*kR);
R = zeros(size(R00));
X = R;
Ntot = maxmodes^2;
pg = 0;
t00 = tic;
for m=1:maxmodes
    for n=1:maxmodes
        t01 = tic;
        if (n==1)&&(m==1)
            % fundamental mode: use analytical solution
            Z = R00+1i*X00;
            pg = pg + 1;
        elseif n>=m
            pg = pg + 1;
            muMax = max(bz(n),bz(m));            
            if useHFapprox 
                intId = find(kR < (muMax*2.5+1));
                intHf = find(kR >= (muMax*2.5));
                if (n==m)
                    R(intHf) = (R00(intHf) .* k(intHf) ./ sqrt(k(intHf).^2 - (muMax./a).^2));
                else
                    R(intHf) = ((R00(intHf)-1) ./ (1 - (muMax./kR(intHf)).^2));
                end
                X(intHf) = (X00(intHf) ./ (1 - (muMax./kR(intHf)).^2));
            else
                intId = 1:length(kR);
            end
            R(intId) = integrateR(n,m,kR(intId),bz,tol);
            X(intId) = integrateX(n,m,kR(intId),bz,1e-6);
            Z = R+1i*X;
        else
            Z = Zmat(n,m,:);
            pg = pg + 1;
        end
        if (progressReport && ~(n<m))
            fprintf('%5.1f%%: Calculated mode (%d,%d) (%.1fs, %.1fs tot)\n',100*pg/Ntot,n,m,toc(t01),toc(t00));
        end
        Zmat(m,n,:) = Z;
    end
end

Zmat = rho*c/S*Zmat;

% Looks like these integration functions can't be vectorized, but if
% someone knows how to, please let me know! The integration is quite a bit
% slower in Octave. 
function R = integrateR(n,m,kR,bz,tol)
R = zeros(size(kR));
if isOctave()
    for nk = 1:length(kR)
        R(nk) = quadgk(@(x) MPM_ASbresistance(x,n,m,kR(nk),bz), 0,pi/2, tol*100);
    end        
else
    for nk = 1:length(kR)
        R(nk) = quadgk(@(x) MPM_ASbresistance(x,n,m,kR(nk),bz), 0,pi/2,'AbsTol', tol);
    end
end

function X = integrateX(n,m,kR,bz,tol)
X = zeros(size(kR));
if isOctave()
    for nk = 1:length(kR)
        X(nk) = quadgk(@(x) MPM_ASbreactance(x,n,m,kR(nk),bz), 0,10, tol);
    end
else
    for nk = 1:length(kR)
        X(nk) = quadgk(@(x) MPM_ASbreactance(x,n,m,kR(nk),bz), 0,10, 'AbsTol', tol);
    end
end

function dR = MPM_ASbresistance(phi,n,m,kR,bz)
sinphi = sin(phi);
dR = sinphi.*FuncDn(sinphi, bz(n), kR).*FuncDn(sinphi, bz(m), kR);

function dX = MPM_ASbreactance(phi,n,m,kR,bz)
coshphi = cosh(phi);
dX = coshphi.*FuncDn(coshphi, bz(n), kR).*FuncDn(coshphi, bz(m), kR);

function D = FuncDn(tau, gamma_n, kR)
D = -sqrt(2).*tau.*besselj(1,(tau.*kR))./ ((gamma_n./kR).^2-tau.^2);


