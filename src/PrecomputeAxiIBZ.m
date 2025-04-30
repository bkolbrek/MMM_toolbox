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

function PrecomputeAxiIBZ(maxModes)
load MMM_besselzeros.mat;
if nargin < 1
    maxModes = 32;
end
fprintf('Precomputing radiation impedance using %d modes\n', maxModes);
disp('This may take several minutes, depending on your computer and the number of modes desired.');
kamax = bz(maxModes)*2.5;
kamin = 0.1;
nk = bz(maxModes) * 4;
filename = sprintf('ZradAS%d.mat', maxModes);

ka1 = logspace(log10(kamin), log10(3), 50);
ka2 = linspace(3, kamax, nk);
ka = [ka1 ka2(2:end)];
% ka = linspace(kamin, kamax, nk);
k = ka*sqrt(pi);
tic;
Zmat = MMM_ASbaffledradzmatrix(k, 1, 1, 1, maxModes, bz, true, true);
toc;
if isOctave()
    str = sprintf('save "-mat" "ZradAS%d.mat" "ka" "Zmat"', maxModes);
    eval(str);
    save "-mat" "zradfile.mat" "filename";
else
    save(filename, 'ka','Zmat');
    save('zradfile.mat', 'filename');
end