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
% Version  : Alpha RA 1.0
%
%########################################################################
%
%	bandpass_filter.m:
%
% Description:
%
%   A  MATLAB class to bandpass filter signals
%

classdef bandpass_filter
    properties
        fs_;
        low_;
        high_;
        low_a_;
        low_b_;
        high_a_;
        high_b_;
    end
   
   methods
       
       function obj = bandpass_filter(fs,low,high)
           
          if nargin > 0
             if isnumeric(fs)
                obj.fs_ = fs;
                obj.low_ = high;
                obj.high_ = low;
                
                fNorm = obj.low_ / (obj.fs_/2);                    
                [obj.low_b_,obj.low_a_] = butter(3, fNorm, 'low'); 
                
                fNorm = obj.high_ / (obj.fs_/2);                    
                [obj.high_b_,obj.high_a_] = butter(3, fNorm, 'high'); 
             else
                error('fs must be numeric')
             end
          end
          
       end
       
       function signal_f = get(obj,signal)
           signal_f = filtfilt(obj.low_b_, obj.low_a_, signal);
           signal_f = filtfilt(obj.high_b_, obj.high_a_, signal_f); 
       end
   end
end

