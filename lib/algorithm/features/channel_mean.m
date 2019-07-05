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
%	skin.m:
%
% Description:
%
%   The classical channel expectation
%

classdef channel_mean
    
   methods
       
       function obj = channel_mean()
       
       end
       
       function [mu, obj] = get(obj,rgb)
           mu=mean(rgb);
       end
   end 
end

