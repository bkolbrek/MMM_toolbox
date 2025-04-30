%MMM_ASradiatedPressure : Calculate the radiated pressure
%
% data = MMM_ASradiatedPressure(data, fieldPoints, useFarfieldApprox)
%
% Calculates the pressure radiated from a horn or other circular radiator where the
% surface volume velocity is described by modes.
%
% Input parameters:
%  data : MMM parameter struct.
%  field points : 
%  useFarfieldApprox : (boolean, optional) selects the type of method used
%       to calculate the pressure. Setting this to true uses a fast
%       far-field approximation
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
function data = MMM_ASradiatedPressure(data, fieldPoints, useFarfieldApprox)

if nargin < 3
    useFarfieldApprox = true;
end

if useFarfieldApprox
    data.pRad = modalRadiatedPressure(data, fieldPoints);
else
    data.pRad = rayleighRadiatedPressure(data, fieldPoints);
end


function prext = modalRadiatedPressure(data, fieldPoints)
a = sqrt(data.Sm/pi);
prext = zeros(size(fieldPoints,1), data.nfreq);
for ii = 1:size(fieldPoints,1)
    pe = fieldPoints(ii,:);
    R = norm(pe);
    theta = atan2(pe(1),pe(2));
    s = data.k*a*sin(theta);
    sm = s(ones(data.nModes,1),:);
    bzm = data.eigenValues(:,ones(1,length(s)));
    Theta2M = 2*sm.*besselj(1,sm)./(sm.^2 - bzm.^2);
    if theta == 0
        Theta2M(1,:) = ones(1,length(s));
    end
    if data.nModes > 1
        ModalSum = sum(1i*(Theta2M.* data.Umouth))';
    else
        ModalSum = 1i*Theta2M.';
    end
    pf = data.rho*data.c/(2*pi*R)*exp(-1i*data.k'*R).*data.k';
    prext(ii,:) = pf.*ModalSum;
end

function prext = rayleighRadiatedPressure(data, fieldPoints)
a = sqrt(data.Sm/pi);
r = (0:a/(data.nIntegrationPoints-1):a);
rp2=r(2:end);
rp1=r(1:end-1);
rp=(rp1+rp2)/2;
phi = MMM_ASgeteigenfunctions(a, rp,  data.eigenValues, true);
Ur = phi*data.Umouth; %volume velocity as function of radius
uo = Ur / (a^2*pi); % mouth particle velocity

prext = zeros(size(fieldPoints,1), data.nfreq);
Np = size(fieldPoints);
for ii = 1:Np(1)
    point = fieldPoints(ii,:);
    pr = MMM_ASrayleighint(data.k, r, uo, point, data.rho, data.c);
    prext(ii,:) = pr;
end


function prext = MMM_ASrayleighint(k, vert, vvel, pe, rho, c)

phiext = 0;
p = [pe(1) 0 pe(2)]';
NV = length(vert);
S = 0; 
nk = length(k);

QA = vert(1:NV-1);
QB = vert(2:NV);
radmid = (QA+QB)/2;
glen = abs(QA-QB);
cirmid = 2*pi*radmid;
NT = ceil(1+cirmid./glen);


for ir=1:NV-1
    NTi = NT(ir);
    theta = (0:NTi-1)*2*pi/NTi;
    dtheta = 2*pi/NTi;
    dS = 0.5*dtheta*(QB(ir)^2-QA(ir)^2);
    q = zeros(3,NTi);
    q(1,:) = radmid(ir)*cos(theta);
    q(2,:) = radmid(ir)*sin(theta);
    q(3,:) = 0;
    S = S+dS*NTi;
    r = q-p;
    r = sqrt(r(1,:).^2+r(2,:).^2+r(3,:).^2)';
    rm = repmat(r, 1, nk);
    phiext = phiext + sum((vvel(ir,:).*dS.*exp(-1i*r*k) ./ rm));
end

prext = 1i*k*rho*c./(2*pi).*phiext;