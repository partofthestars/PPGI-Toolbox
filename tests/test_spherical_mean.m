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
%	test_spherical_mean.m:
%
% Description:
%
%   test of the spherical mean feature extraction on given sample rgb data
%
% 
clear all;
close all;

load('./../media/data/example_data.mat');

if ~exist('skin_pixels')
    disp('error: no skin pixels available. execute test_skin.m first!');
    return; 
end

spm=spherical_mean();

for f=1:size(skin_pixels,2)
    f
    [signal(f,:) ssr]=spm.get(skin_pixels{f});
end

fs=25;
low_frequency=0.5;
high_frequency=2.5;
bpf=bandpass_filter(fs,low_frequency,high_frequency);
signal_filtered=bpf.get(signal);

[pearson, rmse, snr, snr_var, bpm] = ground_truth_stats.get(ppg,signal_filtered(:,1),fs);
