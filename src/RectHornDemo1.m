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
%  HOW TO SET THE NUMBER OF SEGMENTS AND MODES
%   The best way to this is by checking for convergence.
%
%   If there are too few segments, the throat impedance will show wide
%   variations at high frequencies, while for most horns the normalized
%   impedance should flatten out towards unity. If the throat impedance
%   flattens out and then suddenly go crazy, this is because numerical
%   errors are introduced because the segments are too long.
%
%   The required number of modes can be seen from the polar map. The
%   combination of high freqencies and large angles is the most
%   challenging, therefore the upper right hand corner of the polar map
%   will converge the slowest. Before convergence there will be null lines
%   comparable to that of a rigid piston/plane wave. Run the demo with N=1
%   to see what this looks like. As the number of modes is increased, the
%   piston pattern will be pushed towards the upper right hand corner and
%   eventually disappear when the number of modes is sufficient for the
%   horn geometry and frequency range. There may still be some changes, so
%   check for convergence to be safe.
%
%		For rectangular horns, it may be difficult to get enough modes for
%		convergence, because the size of the matrices are much larger, since
%		modes in both directions need to be taken into account. This will slow
%		down calculations.
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
% horn parameters
sth = 10e-4; % throat area
smh = 1500e-4; % mouth area;
Lh = 30e-2; % horn length
aspectThroat = 1.0;
aspectMouth = 1.3; % mouth aspect ratio
dz = Lh/250; % segment length
Tn = 7; % horn parameter (hypex/bessel)
HornTypeH = 'flared conical';
HornTypeV = 'exponential';
AddRadius = false; % set to true to gat a radiused flaring of the mouth
RadR = 0.1; % radius of mouth radius
RFta = 80; % end tangent angle of radius (degrees)

% Create horn contour using the axisymmetric profile generator
xt = sqrt(sth*aspectThroat);
yt = sqrt(sth/aspectThroat);
xm = sqrt(smh*aspectMouth);
ym = sqrt(smh/aspectMouth);

coords1 = MMM_1Dhorncoord(HornTypeH, xt, xm, Lh, Tn, dz, AddRadius, 0, RadR, RFta);
coords2 = MMM_1Dhorncoord(HornTypeV, yt, ym, Lh, Tn, dz, AddRadius, 0, RadR, RFta);
horncoords = [coords1 coords2(:,2)];

figure(1)
plot(horncoords(:,1), horncoords(:,2), horncoords(:,1), horncoords(:,3));


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
N = 25; %maximum number of modes (total, xdirection x ydirection)
c = 344; %sound speed
rho = 1.205;% air density
freq = logspace(log10(fmin), log10(fmax), Nf);

% --- INITIALIZATION
% Must be called every time the above parameters change
data = MMM_init(freq, N, horncoords, 'rect', rho, c);



% plot horn contour
MMM_REplotHorn(data,1)


%% -----------------------------------------------
%  The performance calculations

t00 = tic;
disp(['Calculating radiation impedance']);
% NOTE: this is a hack until I get the proper Zrad calculation implemented.
a = data.steppedCoords(end,2); 
b = data.steppedCoords(end,3);
n = ceil(sqrt(data.nModes));
data.Zrad = MMM_RECbaffledradzmatrixIntp(data.k, data.rho, data.c, a, b, n, data.modeIndex);
data.Zrad = data.Zrad(1:data.nModes, 1:data.nModes, 1:data.nfreq);
toc(t00);

disp(['Calculating horn matrices']);
data = MMM_calculateMatrices(data, false);
toc(t00);
MMM_PlotZth(data,2);


%%
% Default excitation is unit throat velocity: data.Umouth = data.UmouthPw*data.St;
% To change exitation, data.Umouth must be calculated from data.Umouth using a different scale factor.

fieldPoints = [Rext*sin(Angext/180*pi) zeros(size(Angext)) Rext*cos(Angext/180*pi);
	zeros(size(Angext)) Rext*sin(Angext/180*pi) Rext*cos(Angext/180*pi)];

disp('Calculating field point pressure');
tic
data = MMM_REmodalradiatedpressure(data, fieldPoints);
toc
pext = data.pRad;
toc(t00);

nfp = length(Angext);
ia = find(mod(Angext, 10)==0);

figure(3);
subplot(2,1,1)
semilogx(data.fvec, 94+20*log10(abs(data.pRad(ia,:))));
xlim(data.fvec([1,end]));
ylabel('dB SPL');
xlabel('Hz');
title('Field point pressures, horizontal');
grid;

subplot(2,1,2)
semilogx(data.fvec, 94+20*log10(abs(data.pRad(nfp+ia,:))));
xlim(data.fvec([1,end]));
ylabel('dB SPL');
xlabel('Hz');
title('Field point pressures, vertical');
grid;

Lp = 20*log10(abs(data.pRad));
Lp = Lp - Lp(1,:);

figure(4)
subplot(2,1,1);
contourf(data.fvec, Angext, Lp(1:nfp,:), 15);
set(gca, 'xscale', 'log');
xlim(data.fvec([1,end]));
ylabel('Degrees');
xlabel('Hz');
t = 'Polar map, horizontal';
title(t);
colorbar;

subplot(2,1,2);
contourf(data.fvec, Angext, Lp((1+nfp):end,:), 15);
set(gca, 'xscale', 'log');
xlim(data.fvec([1,end]));
ylabel('Degrees');
xlabel('Hz');
t = 'Polar map, vertical';
title(t);
colorbar;




