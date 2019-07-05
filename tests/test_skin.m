%########################################################################
%
%   - PPGI Toolbox - 
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
%	test_skin.m:
%
% Description:
%
%   test of the skin segmentation on given sample rgb data
%

clear all;
close all;

load('./../media/data/example_data.mat');

frames=size(rgb,2);

for f=1:frames
    f
    skin_pixels{f}=skin.get(rgb{f});
end

save('./../media/data/example_data.mat','fs','ppg','rgb','skin_pixels','-v7.3');
