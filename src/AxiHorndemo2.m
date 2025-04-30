%MPM_Horndemo2
%
%  This file contains sample code showing how the MMM toolbox can be used to
%  calculate pressure field inside and outside a horn at a single frequency.
%  By adjusting the number of modes, one can examine where and how the sound field
%  changes, as it approaches convergence. If too few modes are specified,
%  the plot will not represent the actual sound field inside the horn. As
%  the number of modes is increased, the field converges, and the
%  equi-pressure contours become less wiggly.
%  A wiggly equi-pressure contour does not necessarily imply low accuracy
%  in the far field pressure.
%
%  Several horn types can be specified, see MMM_AShorncoord for details.
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

clearvars;
% load the precalculated zeros of Bessel function J1
load 'MMM_besselzeros.mat';

% horn parameters
sth = 10e-4; % throat area
sm = 1500e-4; % mouth area;
Lh = 49e-2; % horn length
dz = Lh/250; % segment length
Tn = 1; % horn parameter (hypex/bessel)
HornType = 'tractrix';
AddRadius = false; % set to true to gat a radiused flaring of the mouth
RadR = 0.1; % radius of mouth radius
RFta = 80; % end tangent angle of radius (degrees)

% Create horn contour
rth = sqrt(sth/pi);
rm = sqrt(sm/pi);
horncoords = MMM_1Dhorncoord(HornType, rth, rm, Lh, Tn, dz, AddRadius, 0, RadR, RFta);

% --- field point parameters
% field point distance
Rext = 3;
% field point angles (degrees)
Angext = linspace(0,90,181)';
% use far field approximation - much faster
useFarfieldApprox = true;


% --- simulation parameters
freq = 1500; % plotting frequency
N = 24; %maximum number of modes
c = 344; %sound speed
rho = 1.205;% air density
addNearfield = true;

% --- INITIALIZATION
% Must be called every time the above parameters change
data = MMM_init(freq, N, horncoords, 'axi', rho, c);


% plot horn contour
MMM_ASplotHorn(data,1)

% Calculate and plot sound field
tic;
MMM_ASpressureDistribution(freq, data, addNearfield, 2);
toc;

