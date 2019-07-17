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
%	test_pos.m:
%
% Description:
%
%   test of the the POS algorithm from TU Eindhoven
%

clear all;
close all;

load('./../media/data/example_data.mat');

if ~exist('skin_pixels')
    disp('error: no skin pixels available. execute test_skin.m first!');
    return; 
end

pos=projection_orthogonal_to_skin(45);
pulse=pos.get(skin_pixels);

fs=25;
low=0.5;
high=2.5;
bpf=bandpass_filter(fs,low,high);
pulse_f = bpf.get(pulse);


[pearson, rmse, snr, snr_var, bpm] = ground_truth_stats.get(ppg,pulse_f,fs);

