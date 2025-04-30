%MMM_ASbaffledradzmatrix : Modal radiation impedance matrix
%
% Zmat = MMM_ASbaffledradzmatrix(k, rho, c, S, maxmodes, bz)
%
% Calculates the modal radiation impedance matrix for the
% end of a circular tube terminated in an infinite baffle
% by interpolation of a lookup table.
%
% The matrix Zmat is a square, symmetrical matrix.
%
% Input parameters:
% k : wavenumber
% rho : density of medium
% c : sound speed in medium
% S : area of tube opening
% maxmodes : maximum number of modes calculated
% filename : (optional): the filename of the precomputed data. If not
%  supplied, the filename stored when precalculating the data will be used.
%  
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

function ZmatOut = MMM_ASbaffledradzmatrixIntp(k, rho, c, S, maxModes, filename)
if nargin < 6
    try
        load('zradfile.mat','filename');
    catch
        error('Error in MMM_ASbaffledradzmatrixIntp: Could not load file "zradfile.m".\nPlease make sure you have precalculated the radiation impedance using MMM_ASbaffledradzmatrixPrecompute.m before calling this function.');
    end
end
try
    load(filename, 'ka', 'Zmat');
catch
    error('Error in MMM_ASbaffledradzmatrixIntp: Could not load file "%s"\nPlease make sure you have precalculated the radiation impedance using MMM_ASbaffledradzmatrixPrecompute.m before calling this function.',filename);
end
if (maxModes > size(Zmat,1))
    error('Higher number of modes requested than precalculated in %s',filename);
end
maxModes = min(size(Zmat,1),maxModes);
ZmatOut = zeros(maxModes, maxModes, length(k));
load MMM_besselzeros.mat;
a = sqrt(S/pi);
kain = k*a;
R00 = 1-besselj(1,2*kain)./(kain);
X00 = 2*struveH1(2*kain)./(2*kain);
intpId = find((kain >= min(ka)) & (kain <= max(ka)));
posId = find(kain > max(ka));
for m=1:maxModes
    for n=1:maxModes
        if (m==1) && (n==1)
            Zmn = R00 + 1i*X00;
        elseif (m <= n)
            minka = 1;
            intpIdR = find((kain >= minka) & (kain <= max(ka)));            
            bzq = max(bz(m), bz(n));
            Y = squeeze(real(Zmat(m,n,:)));
            zp = getPolyCoeff(m,n,bz);
            R1 = polyval(zp, kain(kain < minka).^2)';
            R2 = interp1(ka(:), Y, kain(intpIdR)','spline');
            if (m==n)
                R3 = (R00(posId) .* k(posId) ./ sqrt(k(posId).^2 - (bz(m)./a).^2)).';
            else
                R3 = ((R00(posId)-1) ./ (1 - (bzq./kain(posId)).^2)).';
            end
            %R3 = interp1(ka', Y, kain(posexpol)','linear','extrap');
            R = [R1; R2; R3];
            Y = squeeze(imag(Zmat(m,n,:)));
            X1 = Y(1)*kain(kain < ka(1))'/ka(1);
            X2 = interp1(ka(:), Y, kain(intpId)','spline','extrap');
            X3 = (X00(posId) ./ (1 - (bzq./kain(posId)).^2)).';
            X = [X1; X2; X3];
            Zmn = R +1i*X;
        else
            Zmn = ZmatOut(n,m,:);
        end     
        ZmatOut(m,n,:) = Zmn;
        
    end
end
ZmatOut = rho*c/S*ZmatOut;

function p = getPolyCoeff(n,m,bz)
nlow = min(n,m);
nhigh = max(n,m);
if (nlow == 1) && (nhigh > 1)
    p = zeros(1,4);
    p(2) = -1/(3*bz(nhigh)^2);
    p(1) = -4/(15*bz(nhigh)^4) + 1/(15*bz(nhigh)^2);
else
    p = zeros(1,5);
    p(2) = 4/(15*bz(m)^2*bz(n)^2);
    p(1) = 8/(35*bz(m)^4*bz(n)^2) + 8/(35*bz(m)^2*bz(n)^4) - 2/(35*bz(m)^2*bz(n)^2);
end
    


