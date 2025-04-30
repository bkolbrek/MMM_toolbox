# MMM_toolbox

Repository for the Mode Matching Method (MMM) Toolbox by Bjørn Kolbrek.

Copyright (C) 2012-2025 by Bjørn Kolbrek

## Purpose

The MMM Toolbox simulates horns, either loudspeaker horns or musical instruments, using the Mode Matching Method (MMM)
which describe the soundfield inside a duct by a sum of basis functions (modes). 

The code in this toolbox was written during my PhD and has since been somewhat modified to make it easier to use.

Currently only axisymmetric horns can be simulated, but code for rectangular horns is underway. 

## How to set the number of segment and modes

The best way to this is by checking for convergence. 

If there are too few segments, the throat impedance will show wide
variations at high frequencies, while for most horns the normalized
impedance should flatten out towards unity. If the throat impedance
flattens out and then suddenly go crazy, this is because numerical
errors are introduced because the segments are too long. 

The required number of modes can be seen from the polar map. The
combination of high freqencies and large angles is the most
challenging, therefore the upper right hand corner of the polar map
will converge the slowest. Before convergence there will be null lines
comparable to that of a rigid piston/plane wave. Run the demo with N=1
to see what this looks like. As the number of modes is increased, the
piston pattern will be pushed towards the upper right hand corner and
eventually disappear when the number of modes is sufficient for the
horn geometry and frequency range. There may still be some changes, so
check for convergence to be safe.

Looking at the soundfield inside the horn is another option. 
By adjusting the number of modes, one can examine where and how the sound field
changes, as it approaches convergence. If too few modes are specified,
the plot will not represent the actual sound field inside the horn. As
the number of modes is increased, the field converges, and the
equi-pressure contours become less wiggly.
A wiggly equi-pressure contour does not necessarily imply low accuracy
in the far field pressure.


## License

The MPM Toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 2 of the License, or (at your option) any later version.

The MMM Toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the
MMM Toolbox. If not, see <http://www.gnu.org/licenses/>.
