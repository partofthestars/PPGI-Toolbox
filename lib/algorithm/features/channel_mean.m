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
%	skin.m:
%
% Description:
%
%   The classical green channel spatial expectation of pixel intensities
%
% References:
%
%	M. Hülsbusch. A functional imaging technique for opto-electronic assessment of skin perfusion. 
%	PhD thesis, RWTH Aachen University, 2008.
%
%	Verkruysse, W., Svaasand, L. O., & Nelson, J. S. (2008). Remote plethysmographic imaging using ambient light. 
%	Optics express, 16(26), 21434-21445
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

