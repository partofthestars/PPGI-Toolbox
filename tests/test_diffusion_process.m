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
%	test_diffusion_process.m:
%
% Description:
%
%   test of the Diffusion process model. 
%

clear all;
close all;

load('./../media/data/example_data.mat');

if ~exist('skin_pixels')
    disp('error: no skin pixels available. execute test_skin.m first!');
    return; 
end

for f=1:size(skin_pixels,2)
    f
    raw_mean(f,:)=mean(skin_pixels{f});
end

%pre-filter to make it a bit easier for the state space model
fs=25;
low=0.5;
high=2.5;
bpf=bandpass_filter(fs,low,high);
raw_mean_f = bpf.get(raw_mean);

%diffusion process (btw. dynamic bayesian state space model (dbssm) )
%

%estimate of measurement noise standard deviation
%in the IMM model.
sd = 0.025;
%the process spectral density for the bias model rep-
%resents the continuous time noise in the sensor signal
bq = 0.01;
%the resonator process noise spectral density defines
%the continuous-time variation of the resonator sig-
%nals. adjust primarily this parameter to control the
%behavior of the periodic signals.
qr = 0.000001;

%number of harmonics including fundamental
nharm=1;
%time delta
dt=1/fs;
%transition probability between consecutive steps of
%frequencies (i.e. the probability of a jump from e.g. 70 bpm to 71 bpm)
ptrans = 0.0015;
%transition probability between all steps of frequencies
%(i.e. the probability of a jump from e.g. 70 bpm to 80 bpm).
poverall = 0;
%frequency search space
freqlist=45:120;

dbssm=diffusion_process();

pulse=dbssm.get(raw_mean_f(:,2).*0.001,dt,freqlist,nharm,bq,sd,qr,ptrans,poverall);

[pearson, rmse, snr, snr_var, bpm] = ground_truth_stats.get(ppg,pulse,fs);

figure;
plot(pulse*10000);
hold on;
plot(raw_mean_f(:,2),'r')

