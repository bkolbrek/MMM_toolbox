function [plotcoordszOut, plotcoordsxOut, PmatxOut]...
    = MMM_ASpressureDistribution(freq, data, addNearfield, plotOn, resolution)
% 
% Plots the sound field inside and in fromt of a horn. 
% 
% Input parameters:
%  k : single frequency wave number
%  data : MMM data struct
%  addNearfield : adds a region in front of the horn where the nearfield
%       radiated pressure is calculated.
%  ploton : plot the pressure distribution if larger than zero, also the
%   figure number. Uses new figure as default.
%  resolution : the number of points radially, effectively the density of
%       the mesh. Default: 30 points.
%
% Output parameters:
%  plotcoordszOut, plotcoordsxOut : mesh coordinates
%  PmatxOut : complex pressure at the mesh points
%
% See plotting function below for an example of how to use this data.
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
if nargin < 5
    resolution = 30;
end

if nargin < 3
    addNearfield = false;
end

U0 = zeros(data.nModes, 1);
U0(1) = data.St;

data.keepZmatrix = true;
data.fvec = freq;
data.nfreq = 1;
k = freq*2*pi/data.c;
data.k = k;

data.Zrad = MMM_ASbaffledradzmatrix(k, data.rho, data.c, data.Sm, data.nModes, data.eigenValues, true);
data = MMM_calculateMatrices(data, false);

Nz = length(data.steppedCoords);
NP = ceil(Nz/2+1);

if addNearfield
    data.nIntegrationPoints = max(30, resolution + 5);
    dz = data.c/(freq*6);
    boxSize = 2*data.steppedCoords(end,2);
    Nfield = max(10, ceil(boxSize/dz));
    mindist = 0;
    NPcoords = NP + Nfield; 
    z = linspace(mindist, mindist+boxSize, Nfield);
    x = linspace(0, boxSize, resolution);
    [Zm, X] = meshgrid(z, x);
%     X(:,1) = linspace(0,data.steppedCoords(end,2),Npoints);
    fieldPoints = [reshape(X,[],1), reshape(Zm,[],1)];
else
    NPcoords = NP;
end
Pmat = zeros(data.nModes, NP);
Pmatx = zeros(resolution,NPcoords);
plotcoordsz = zeros(resolution,NPcoords);
plotcoordsz(:,1:NP) = ones(resolution,1)*data.steppedCoords([1:2:Nz,Nz],1)';
plotcoordsx = zeros(resolution,NPcoords);
if addNearfield
    plotcoordsx(:,NP+1:end) = X;
    plotcoordsz(:,NP+1:end) = Zm+data.steppedCoords(end,1);
end

%% Calculate the horn pressure field
ip = 1;
plotcoordsx(:,ip) = linspace(0,data.steppedCoords(1,2),resolution);
phix = MMM_ASgeteigenfunctions(data.steppedCoords(1,2), plotcoordsx(:,ip)',  data.eigenValues, true);
Pmat(:,ip) = data.BigZ(:,:,1) * U0;
Pmatx(:,ip) = phix*Pmat(:,ip);
for ik = 1:data.nfreq
    U = U0;

    % propagate back to throat
    for iz = 1:(Nz-1)
        R1 = data.steppedCoords(iz,2);
        R2 = data.steppedCoords(iz+1,2);
        L = data.steppedCoords(iz+1,1) - data.steppedCoords(iz,1);
        if (L>0)
            ip = ip + 1;
            krc = data.k(ik)*data.rho*data.c;
            Z = data.BigZ(:,:,iz,ik);
            kn = MMM_ASmakekm(data.k(ik),data.steppedCoords(iz,:),data.nModes,data.eigenValues);                
            D2 = (1i*sin(L*kn));
            Zc = (krc ./ (data.S(iz)*kn));
            invZc = (data.S(iz)*kn)./(krc); 
            E = diag(exp(-1i*L*kn)); 
            U = (-diag(D2.*invZc)*(Z-diag(Zc))+E) * U;
            Pmat(:,ip) = data.BigZ(:,:,iz+1,1)*U;
            plotcoordsx(:,ip) = linspace(0,data.steppedCoords(iz+1,2),resolution);
            phix = MMM_ASgeteigenfunctions(data.steppedCoords(iz+1,2), plotcoordsx(:,ip)',  data.eigenValues, true);
            Pmatx(:,ip) = phix*Pmat(:,ip);
        else
            F = data.bigF(:,:,iz);
            if R1>R2
                U = (F.')\U;
            else
                U = F.'*U;
            end
        end
    end   
end

%% Add the nearfield pressure
if addNearfield
    useFarfieldApprox = false;
    data = MMM_ASradiatedPressure(data, fieldPoints, useFarfieldApprox);
    pRad = reshape(data.pRad, resolution, Nfield);
    Pmatx(:,NP+1:end) = pRad;
    ind = find(plotcoordsx(:,NP+1) <= data.steppedCoords(end,2));
    coords = plotcoordsx(ind,NP+1);
    phix = MMM_ASgeteigenfunctions(data.steppedCoords(iz+1,2), coords',  data.eigenValues, true);
    p = phix*Pmat(:,ip);
    Pmatx(ind,NP+1) = p;
end

if nargout > 0
    plotcoordszOut = plotcoordsz;
end
if nargout > 1
    plotcoordsxOut = plotcoordsx;
end
if nargout > 2
    PmatxOut = Pmatx;
end

%% Plot the results
if nargin < 4
    figure();
elseif plotOn <= 0
    return
else
    figure(plotOn);
    contourf(plotcoordsz', plotcoordsx', 94+20*log10(abs(Pmatx))', 25);
    colorbar;
    axis equal;
    hold on;
    plot([data.rawCoords(:,1); data.rawCoords(end,1)],[data.rawCoords(:,2); 0],'k');
    hold off; 
    title(sprintf('Sound field at %.1fHz, using %d modes', freq, data.nModes));
    xlabel('z [m]');
    ylabel('Radius [m]');
end
%%
    






