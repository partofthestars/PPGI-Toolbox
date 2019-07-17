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
%	test_ssr.m:
%
% Description:
%
%   the spatial subspace rotation algorithm from TU Eindhoven
%

clear all;
close all;

addpath('./../lib/utils')
addpath('./../lib/algos')

load('./../media/data/example_data.mat');

if ~exist('skin_pixels')
    disp('error: no skin pixels available. execute test_skin.m first!');
    return; 
end

fs=25;
window_size=3;
overlap=2;

ssr=spatial_subspace_rotation(fs,window_size,overlap);

frames=size(rgb,2);

for f=1:frames
    f
    [pulse(f,:) ssr]=ssr.get(skin_pixels{f});
end

low=0.5;
high=2.5;
bpf=bandpass_filter(fs,low,high);
pulse_f = bpf.get(pulse);

[pearson, rmse, snr, snr_var, bpm] = ground_truth_stats.get(ppg,pulse_f,fs);