%MMM_1Dhorncoord(hType, Sth, Sm, L, Tn, dz, addradius, thFta, RadR, RadFta)
%
% horncoord = MMM_AShorncoord(hType, Sth, Sm, L, Tn, dz, addradius, thFta, RadR, RadFta)
%
% Calculates a horn contour of a given type, from parameters. An optional
% flared/radiused mouth can be added by setting the addradius parameter to
% true, and supplying values for the optional parameters thFta, RadR and
% RadFta.
%
% The horn is axisymmetric. The coordinates are in (z,y) format.
%
% Input parameters:
% hType : horn type, string, case insensitive. The following horn types are supported:
%   Conical
%   Exponential (or Expo)
%   Hypex (Hyperbolic-exponential / Salmon type horn)
%   OSWG (Oblate Spheroidal Waveguide)
%   Bessel (Bessel horn, S = S0 * x^n)
%   Spherical (Spherical wave horn. To specify cutoff wavenumber, set Sm = -kc)
%   Tractrix (Tractrix horn. To specify cutoff wavenumber, set L = -kc)
%   Radius : (the horn profile is part of a circle or radius RadR, starting at an angle
%       thFta, and ending at an angle RadFta (Fta = flare tangent angle))
%	Flared conical : keele type
% Yth : throat dimension, m
% Ym : mouth dimension, m
% L : horn length, m
% Tn : parameter T for hypex horns, or parameter n for Bessel horns and
%	flared conical
% dz : maximum segment length, m
% addradius : boolean, optional. Adds a radius to the horn.
%    If horntype 7 (radius) is used, this parameter must be set to false, and
%    the remaining optional parameters must be supplied.
%    If a radius should be added to the horn, this parameter must be set to
%    true, and the remaining parameters must be supplied. thFta can be set
%    to an arbitrary value.
% thFta : (optional, see above) throat flare tangent angle, degrees. Start angle for radius type horn.
% RadR : (optional, see above) radius of circle defining radius type horn.
% RadFta : (optional, see above) mouth flare tangent angle, degrees. Max 90 degrees.
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

function horncoord = MMM_1Dhorncoord(hType, Yth, Ym, L, Tn, dz, addradius, thFta, RadR, RadFta)
if nargin < 7
    addradius = false;
end
if L>0
    N = ceil(L/dz);
    dz = L/N;
    horncoord = zeros(N, 2);
end
%Yth = sqrt(Sth/pi);


if (Yth==abs(Ym))
    horncoord(:,1) = linspace(0, L, N);
    horncoord(:,2) = ones(N,1)*Yth;
else
    switch lower(hType)
        case 'conical'
            x1 = L/(Ym/Yth-1);
            horncoord(:,1) = linspace(0, L, N);
            horncoord(:,2) = Yth*((horncoord(:,1)+x1)./x1);
        case {'exponential', 'expo'}
            m = 1/L*log(Ym/Yth);
            horncoord(:,1) = linspace(0, L, N);
            horncoord(:,2) = Yth*exp(horncoord(:,1)*m);
        case 'hypex'
			Sth = Yth^2*pi;
			Sm = Ym^2*pi;
            %m = 1/(2*L)*log(Sm/Sth);
            tt = -2*Sth + 2*Sth*Tn^2 + 4*Sm + 4*sqrt(-Sm*Sth+Sm*Sth*Tn^2+Sm^2);
            tn = 2*Sth*(Tn*(Tn+2)+1);
            m = 1/(2*L)*log(tt/tn);
            horncoord(:,1) = linspace(0, L, N);
            horncoord(:,2) = Yth*(cosh(horncoord(:,1)*m)+Tn*sinh(horncoord(:,1)*m));
        case 'oswg'
            rt2 = Yth^2;
            rm2 = Ym^2;
            tanang2 = (rm2-rt2)/(L^2);
            horncoord(:,1) = linspace(0, L, N);
            horncoord(:,2) = sqrt(rt2+tanang2*(horncoord(:,1).^2));
        case 'bessel'
			Sth = Yth^2*pi;
			Sm = Ym^2*pi;
            x1 = L/((Sm/Sth)^(1/Tn)-1);
            horncoord(:,1) = linspace(0, L, N);
            horncoord(:,2) = Yth*((horncoord(:,1)+x1)./x1).^(Tn/2);
        case 'spherical'
			Sm = Ym^2*pi;
            if (Sm<0) %Sm specifies cutoff wavenumber
                kc = -Sm;
                [~, Xhmax] = CalcSphericalMaxL(Yth, kc);
                if L > Xhmax
                    L = Xhmax;
                    disp(['MMM_1Dhoorncoord: L to large! Value changed to ' num2str(L) 'm']);
                end
                N = ceil(L/dz);
                dz = L/N;
                horncoord = zeros(N, 2);
                horncoord(:,1) = linspace(0, L, N);
            else
                % find cutoff wavenumber
                R1 = Yth;
                R2 = Ym;
                kcmax = log(Ym/Yth)/L;
                kcmin = 0;
                Diff = 1;
                while abs(Diff) > 1e-6
                    kc = (kcmin + kcmax)/2;
                    m = 2*kc;
                    R0 = 4/m;
                    H0 = R0 - sqrt(R0^2 - R1^2);
                    H = R0 - sqrt(R0^2 - R2^2);
                    VX = log(H/H0)/m;
                    if ((H0+VX-H) > L)
                        kcmin=kc;
                    else
                        kcmax=kc;
                    end
                    Diff = kcmax-kcmin;
                end
                
                disp(['Spherical horn cutoff kc = ' num2str(kc) ' ( = ' num2str(kc*344/2/pi) 'Hz)']);
                if abs(imag(VX))>0
                    disp('Invalid horn data - change length!');
                    L = real(H0+VX-H);
                end
                
                N = ceil(L/dz);
                dz = L/N;
                horncoord = zeros(N, 2);
                horncoord(:,1) = linspace(0, L, N);
            end
            for ii=1:N
                horncoord(ii,2) = GetSphericalYx(Yth, kc, horncoord(ii,1));
            end
        case 'tractrix'
            if (Yth>Ym)
                warning('MMM_AShoorncoord: Throat area Sth must be smaller than mouth area Sm for tractrix horn!');
            else
                if L<0 %lenght specifies cutoff wavenumber
                    kc = -L;
                else
                    % find cutoff wavenumber
                    Yx = Ym;
                    kcmax = 1/Yx;
                    kcmin = 0;
                    Xmin = tractrix(Yth,Yx);
                    if L < Xmin
                        % There is a minimum length given by
                        L = Xmin;
                        disp(['MMM_AShoorncoord: L to small! Value changed to ' num2str(L) 'm']);
                        kc = kcmax;
                    else
                        Diff = 1;
                        while abs(Diff) > 1e-6
                            kc = (kcmax+kcmin)/2;
                            Xmax = tractrix(Yth,1/kc);
                            XfromM = tractrix(Yx,1/kc);
                            tx = Xmax - XfromM;
                            Diff = tx - L;
                            if Diff > 0
                                kcmin = kc;
							else
                                kcmax = kc;
                            end
                        end
                    end
                end
                disp(['Tractrix cutoff kc = ' num2str(kc) ' ( = ' num2str(kc*344/2/pi) 'Hz)']);
                Ym = 1/kc;
                Xmax = tractrix(Yth,Ym);
                XfromM = tractrix(Yx,Ym);
                L = real(Xmax-XfromM);
                N = ceil(L/dz);
                dz = L/N;
                horncoord = zeros(N, 2);
                horncoord(:,1) = linspace(0, L, N);
                horncoord(1,2) = Yth;
                for ii = 2:N
                    horncoord(ii,2) = CalcTractrixAtX(Yth, kc, horncoord(ii,1));
                end
            end
        case 'radius'
            RadFta = max(RadFta, 90);
            r0 = Yth - RadR*(1-cos(thFta*pi/180));
            z0 = -RadR * sin(thFta*pi/180);
            zL = RadR * sin(RadFta*pi/180);
            L = z0 + zL;
            N = ceil(L/dz);
            horncoord = zeros(N, 2);
            horncoord(:,1) = linspace(0, L, N);
            theta = real(asin((horncoord(:,1)-z0)/RadR));
            horncoord(:,2) = r0 + RadR - RadR*cos(theta);
		case 'flared conical'
			a = Yth;
			b = (Ym/1.5 - Yth) / L;
			c = (Ym - b*L - a)/(L^Tn);
			x = linspace(0, L, N);
			horncoord(:,1) = x;
			horncoord(:,2) = a+ b*x + c*x.^Tn;
        otherwise
            horncoord(:,1) = linspace(0, L, N);
            horncoord(:,2) = ones(N,1)*Yth;
            warning(['Horn type ' hType ' is not supported.']);
    end
end

% if flag is set, add a radius at the mouth, given by RadR and RadFta
if addradius
    maxs = horncoord(end,2)^2*pi;% max(horncoord(:,2))^2*pi;
    fta1 = 180/pi*atan((horncoord(end,2)-horncoord(end-1,2)) / (horncoord(end,1)-horncoord(end-1,1)));
    horncoord2 = MMM_AShorncoord('radius', maxs, 0, L, Tn, dz, false, fta1, RadR, RadFta);
    horncoord2(:,1) = horncoord2(:,1)+max(horncoord(:,1));
    horncoord = [horncoord ; horncoord2];
end


function y = CalcTractrixAtX(Yt, kc, x)
Ym = 1/kc;
Xmax = tractrix(Yt,Ym);
xfromM = Xmax-x;
Ymax = Ym;
Ymin = Yt;
Diff = 1;
while (abs(Diff) > 1e-8)
    y = (Ymin+Ymax)/2;
    tx = tractrix(y,Ym);
    Diff = tx-xfromM;
    if Diff>0
        Ymin = y;
    else
        Ymax = y;
    end
end

function x = tractrix(Y1, Y2)
x = Y2 * log((Y2 + sqrt(Y2^2 - Y1^2))/Y1) - sqrt(Y2^2 - Y1^2);


function [Xmax, Xhmax] = CalcSphericalMaxL(y0, kc)
m = 2*kc;
r0 = 2/kc;
h0 = r0 - sqrt(r0^2 - y0^2);
Xmax = 1/m * log(1/(m*h0));
Xhmax = Xmax + h0 * (1 - exp(m*Xmax));


function Xvir = GetSphericalVirLength(y0, kc, Lphys)
m = 2*kc;
r0 = 2/kc;
h0 = r0 - sqrt(r0^2 - y0^2);
Xvir_min = Lphys;
Xvir_max = 1/m * log(1/(m*h0));
Xhmax = Xvir_max + h0 * (1 - exp(m*Xvir_max));
if max(Lphys) > Xhmax
    Xvir = Xvir_max;
else
    Diff=1;
    while (abs(Diff) > 1e-8)
        Xvir = (Xvir_min + Xvir_max) / 2;
        Diff = Xvir + h0 * (1 - exp(m * Xvir)) -  Lphys ;
        if Diff < 0
            Xvir_min = Xvir;
        else
            Xvir_max = Xvir;
        end
    end
end

function Yx = GetSphericalYx(y0, kc, x)
Xv = GetSphericalVirLength(y0, kc, x);
Yx = GetSphericalYxVir(y0, kc, Xv);

function Yx = GetSphericalYxVir(y0, kc, xvir)
m = 2*kc;
r0 = 2/kc;
h0 = r0 - sqrt(r0^2 - y0^2);
h = h0 * exp(m*xvir);
Ak = 2*pi*r0*h;
Yx = sqrt(Ak/pi - h^2);

