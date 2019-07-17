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
%	test_bland_altman.m:
%
% Description:
%
%   test of the Blant-Altman and correlation graph 
%   for two datasets.
%

clear all;
close all;

load('./../media/data/example_data.mat');

if ~exist('skin_pixels')
    disp('error: no skin pixels available. execute test_skin.m first!');
    return; 
end

chm=channel_mean();

for f=1:size(skin_pixels,2)
    f
    [signal(f,:) ssr]=chm.get(skin_pixels{f});
end

fs=25;
low_frequency=0.5;
high_frequency=2.5;
bpf=bandpass_filter(fs,low_frequency,high_frequency);
signal_filtered=bpf.get(signal);

[pearson, rmse, snr, snr_var, bpm, bpm_ppg] = ground_truth_stats.get(ppg,signal_filtered(:,2),fs);

% regions per data
territories = {'Green channel mean'};%,'Other Algo','OtherALgo2'};
nterritories = length(territories);

% user states during measurement
states = {'1 users - head movment'};
nstates = length(states);

data1 = bpm_ppg;
data2 = bpm;

% ba plot paramters
tit = 'Heart Rate: PPG vs. PPGI'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'PPGI','PPG','BPM'}; % names of data sets
corrinfo = {'n','r','sse'}; % stats to display of correlation scatter plot
ba_info={''};%ba_info = {'RPC(%)'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 1
	colors = 'br'; % colors for the data sets
else
	colors = [0 0 1;...
		      1 0 0];
end
symbols = ''; % symbols for the data sets (default)

% generate figure with symbols
ba=bland_altman();
[cr, fig, statsStruct]=ba.draw(data1', data2',label,tit,gnames,corrinfo,ba_info,limits,colors,symbols);
