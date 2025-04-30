%MMM_ASgetDI
%
%  DI = MMM_ASgetDI(data, angles, pext)
%
%  Calculates the directivity index from an array of field point pressures.
%  DI is calculated from Beranek, eq. 4.24 and 4.19.
%  If not enough field points are provided, the integration is performed
%  using weights calculated from Gerzon: "Calculating the Directivity Factor
%  of Transducers from Limited Polar Diagram Information", JAES 1975, p.
%  369
%
%  Input parameters:
%  k : wavenumber, assumed to be a row vector
%  angles : angles in degrees, assumed to be a column vector
%  pext : field point pressures, na rows by nk columns
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
function data = MMM_ASgetDI(data, angles)
pMag = abs(data.pRad);
if length(angles)>100
    dtheta = pi/180*angles(2);
    sia = dtheta*sin(pi/180*angles);
    sia = sia(:,ones(1,data.nfreq));
    Wrad = sum((pMag).^2.*sia);
    Q = 2*(pMag(1,:).^2)./Wrad;
else    
    dang = angles(2);
    n = round(180/dang);

    m = n/2 + 1;
    wt = zeros(m,1);
    % calculate weights for numerical integration (after Gerzon 1975)
    for r=0:2:n
        if (r==0)||(r==n)
            k1 = 0.5;
        else
            k1 = 1;
        end
        wt(1) = wt(1) + k1*(-1/(r^2-1));
    end
    wt(1) = wt(1)/(n);

    for i=1:m-1
        wt(i+1) = 0;
        for r=0:2:n
            if (r==0)||(r==n)
                k1 = 0.5;
            else
                k1 = 1;
            end
            wt(i+1) = wt(i+1)+ k1*(-1/(r^2-1)) * cos(pi*r*i/n);
        end
         wt(i+1) = 2*wt(i+1)/n;
    end
    s = size(pMag);
    wt = wt(:,ones(1,s(2)));
    Q = sum(wt.*pMag.^2);
    Q = pMag(1,:).^2./Q;

    maxq = (0.5*n+1)^2*sqrt(2);
    if max(Q) > maxq
        warning('Directivity index larger than %2.0fdB is unreliable. \nPlease use more field points.',10*log10(maxq));
    end  
end
data.DI = 10*log10(abs(Q));