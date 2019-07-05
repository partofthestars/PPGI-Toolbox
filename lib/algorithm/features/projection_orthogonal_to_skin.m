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
%	projection_orthogonal_to_skin.m:
%
% Description:
%
%   Implements the POS algorithm from TU Eindhoven
%
% References:
%
%   Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG.
%   IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491
%

classdef projection_orthogonal_to_skin
    
   properties
       fs_;
       window_size_;
       overlap_;
       shift_;
       U_;
       S_;
       SR_B_;
   end
   
   methods
       
       function obj = projection_orthogonal_to_skin(fs,window_size,overlap)
           
          if nargin > 0
             if isnumeric(fs)
                obj.fs_ = fs;
                obj.window_size_ = window_size;
                obj.overlap_ = overlap;
                obj.shift_ = 1;%obj.window_size_-obj.overlap_;
                obj.U_= circular_buffer(500000,true);
                obj.S_= circular_buffer(500000,true);
                obj.SR_B_= circular_buffer(fs*(window_size+overlap)+500000,false);
             else
                error('fs must be numeric')
             end
          end
       end
       
       function ret =size(obj)
           ret = obj.SR_B_.size();
       end
     
       function [pulse, obj] = get(obj,skin_pixels)
           
           pulse=0;%[0 0 0];

           [rows cols]=size(skin_pixels);

          
       end
       
       %function obj = reset(obj)
       %    obj.U_.reset();
       %end
   end
end

