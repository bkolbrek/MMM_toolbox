%MPM_getimpedances : Calculate the impedances through a horn
%
% BigZ = MPM_getimpedances(data, progressreport)
%
% Calculates the modal impedances at every duct junction in a horn defined
% by the coordinate list coords.
% The F matrices (BigF) must have been calculated on beforehand.
% The resulting matrix BigZ contains the modal impedances n,m at point iz
% and wavenumber index ik as BigZ(n,m,iz,ik)
%
% Input parameters:
% data : MMM parameter struct.
% progressreport : (boolean, optional) prints the current wavenumber and the
%    percentwise progress. Default is off. Calculations are slightly faster
%    with this option turned off.
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
function data = MMM_calculateMatrices(data, progressreport)

if nargin<2
    progressreport = false;
end

nModes = data.nModes;
data.Umat = zeros(data.nModes, data.nModes, data.nfreq);
if (data.keepZmatrix)
    data.BigZ = zeros(nModes,nModes,length(data.steppedCoords), length(data.k));
    data.BigZ(:,:,end,:) = data.Zrad;
end
for ik = 1:data.nfreq
    if progressreport
        disp(['Calculating k = ' num2str(data.k(ik)) ' (' num2str(ik/data.nfreq*100) '%)']);
    end
    U = eye(nModes);
    Z = data.Zrad(:,:,ik);
    % propagate back to throat
    for iz = (size(data.steppedCoords,1)-1):-1:1
        c1 = data.steppedCoords(iz,:);
        c2 = data.steppedCoords(iz+1,:);
        L = data.steppedCoords(iz+1,1) - data.steppedCoords(iz,1);
        if (L>0) %propagate along straight duct
            krc = data.k(ik)*data.rho*data.c;
            kn = data.makekm(data.k(ik), c1, data.nModes, data.modeInfo);
            D2 = (1i*sin(L*kn));
            D3 = (tan(L*kn));
            Zc = (krc./(data.S(iz)*kn));
            D2Zc = diag(Zc./D2); % both matrices are diagonal
            iD3Zc = diag(Zc./(1i*D3));
            Z = iD3Zc - D2Zc/(Z+iD3Zc) * D2Zc;
            D2 = (1i*sin(L*kn));
            invZc = (data.S(iz)*kn)./(krc);
            E = diag(exp(-1i*L*kn));
            U = U*(-diag(D2.*invZc)*(Z-diag(Zc))+E);
        else %propagate across discontinuety
            F = data.bigF(:,:,iz);
            Ft = F.';
            if MMM_useV(c1, c2)
                Z = F\Z/Ft;
                U = U/(Ft);
            else
                Z = F*Z*Ft;
                U = U*Ft;
            end
        end
        if (data.keepZmatrix)
            data.BigZ(:,:,iz,ik) = Z; %keep the impedance for sound field calculation
        end
    end
    data.Umat(:,:,ik) = U;
end

% get throat impedance
data.Z00 = squeeze(data.BigZ(1,1,1,:));

% get mouth volume velocity matrix for throat plane wave
data.UmouthPw = squeeze(data.Umat(:,1,:));
% default volume
data.Umouth = data.UmouthPw*data.St;
