%########################################################################
%
%	- PPGI Toolbox - 
%   A MATLAB toolbox for Photoplethysmography Imaging (PPGI)
%
% Author   : Christian S. Pilz
% Company  : The Nature of Space of Time
% Date     : 07.05.2019
%
% Contact  : cpi@partofthestars.com
% Web Page : www.partofthestars.com
%
% Version  : beta0.1
%
%########################################################################
%
%	test_stochastic_resonator.m:
%
% Description:
%
% 	The frequency trajectories are drawn randomly such that the
% 	random walk (Wiener process) defines a trajectory that is
% 	transformed into a frequency time series f (in Hertz). The drifter
% 	method is applied to the data after the simulated signal has been
% 	transformed to downsampled observations (defined by dts). The mean
% 	squared error results are then captured and visualized.
%

%
close all; %% Simulate M draws

% Number of random draws
M = 1;

% Different time discretizations (TR) to consider
dts = [0.01 0.05 0.1:0.1:1 1.2:.2:2.4];
% Allocate space for results
MSE = zeros(M,numel(dts));
C = zeros(M,numel(dts));

for i=1:M
	[MSE(i,:),C(i,:),x,f] = simulate_and_estimate(dts,25);
	figure;
	subplot(2,1,1)
	plot(x(1,:),'black');
	title('Stochastic Oscillator');
	ylabel('Amplitude');
	xlabel('Time in seconds')
	x_ticks = 0:0.01:25;
	set(gca,'XTickLabel',x_ticks(1:2500).*500 );
	subplot(2,1,2)
	plot(f,'black')
	ylabel('Frequency in Hz');
	xlabel('Time in frames')
	title('Oscillator Frequency')
end
