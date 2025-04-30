%MPM_Horndemo1
%
%  This file contains sample code showing how the MMM toolbox can be used to
%  calculate the performance of a horn exited by a plane wave at the throat.
%  The horn contour is first calculated, and cast into the form of a stepped
%  duct. Next, the modal radiation impedance at the mouth is calculated by 
%  interpolation to speed up the calculations.
%  Following this, the throat impedance, and the modal impedances throughout
%  the horn are calculated. This data is then used in the propagation of
%  velocity from throat to mouth. The mouth velocity is used to compute the
%  resulting pressure in the free field, using a far field modal approach. 
%
%  Several horn types can be specified, see MMM_AShorncoord for details. 
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

clearvars;
% load the precalculated zeros of Bessel function J1

 
% horn parameters
sth = 10e-4; % throat area
sm = 500e-4; % mouth area;
Lh = 30e-2; % horn length
dz = Lh/250; % segment length
Tn = 1; % horn parameter (hypex/bessel)
HornType = 'exponential'; 
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
fmin = 100; %lower frequency
fmax = 15000; % upper frequency
Nf = 200; % number of frequencies
N = 8; %maximum number of modes
c = 344; %sound speed
rho = 1.205;% air density
freq = logspace(log10(fmin), log10(fmax), Nf);

% --- INITIALIZATION 
% Must be called every time the above parameters change
data = MMM_init(freq, N, horncoords, 'axi', rho, c);

 
% plot horn contour
MMM_ASplotHorn(data,1)


%% -----------------------------------------------
%  The performance calculations

t00 = tic;
disp(['Calculating radiation impedance']);
data.Zrad = MMM_ASbaffledradzmatrixIntp(data.k,data.rho,data.c,data.Sm,data.nModes);
toc(t00);

disp(['Calculating horn matrices']);
data = MMM_calculateMatrices(data, false);
toc(t00);
MMM_PlotZth(data,2);

fieldPoints = [Rext*sin(Angext/180*pi) Rext*cos(Angext/180*pi)];

%%
% Default excitation is unit throat velocity: data.Umouth = data.UmouthPw*data.St;
% To change exitation, data.Umouth must be calculated from data.Umouth using a different scale factor.
disp(['Calculating field point pressure']);
tic
data = MMM_ASradiatedPressure(data, fieldPoints, useFarfieldApprox);
toc
pext = data.pRad;
toc(t00);

figure(3);
ia = find(mod(Angext, 10)==0);
semilogx(data.fvec, 94+20*log10(abs(data.pRad(ia,:))));
xlim(data.fvec([1,end]));
ylabel('dB SPL');
xlabel('Hz');
%legend(num2str(Angext(ia)), 'Location', 'NorthWest')
title('Field point pressures');
grid;

MMM_ASpolarMap(data, Angext, true, true, 4);


