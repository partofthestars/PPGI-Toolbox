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
%	local_group_invariance.m:
%
% Description:
%
%   Implements the Local Group Invariance algorithm.
%   NOTE: the implementation considers the translation invariance solely.
%
% References:
%
%	Christian S. Pilz, S. Zaunseder, J. Krajewski, V. Blazek, 
%   Local Group Invariance for Heart Rate Estimation from Face Videos in the Wild, 
%   The IEEE Conference on Computer Vision and Pattern Recognition (CVPR) Workshops, 
%   pp.1254-1262, Salt Lake City, 2018
%

classdef local_group_invariance
    
   properties
 
   end
   
   methods
       
       function obj = local_group_invariance()
          
       end
      
       function [pulse, obj] = get(obj,skin_pixels)
           
           pulse=[];
           frames=size(skin_pixels,2);
           C=[];
           
           for f=1:frames
               x=skin_pixels{f};
               C(f,:)=[mean(x(:,1)); mean(x(:,2)); mean(x(:,3));];
           end
           
           center = mean(C);
           centered = bsxfun(@minus,C,center);
           
           [U,E,V]=svd(centered');
           
           S=U(:,1)';
           P=eye(3)-S'*S;%rank 1
           
           f_T(1:frames,3)=0;
           for f=1:frames
              f_T(f,:)=(P*centered(f,:)')';   
           end

           pulse=double(f_T(:,2));
       end
   end
end

