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
%	projection_orthogonal_to_skin.m:
%
% Description:
%
%   Implements the POS algorithm from TU Eindhoven
%
% References:
%
%	Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG.
%   IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491
%

classdef projection_orthogonal_to_skin
    
   properties
       window_size_;
   end
   
   methods
       
       function obj = projection_orthogonal_to_skin(window_size)
           
          if nargin > 0
             if isnumeric(window_size)
                obj.window_size_ = window_size;
             else
                error('fs must be numeric')
             end
          end
       end
      
     
       function [pulse, obj] = get(obj,skin_pixels)
           
           pulse=[];
           frames=size(skin_pixels,2);
           C=[];
           
           for f=1:frames
               x=skin_pixels{f};
               C(:,f)=[mean(x(:,1)); mean(x(:,2)); mean(x(:,3));];
           end

           L=obj.window_size_;
           H(1,1:frames)=0;
           for f=1:frames-L+1
               block=C(:,f:f+L-1)';
               mu_C=mean(block,1);
               C_normed=[];
               for t=1:size(block,1)
                   C_normed(t,:)=block(t,:)./mu_C;
               end
               S=[0,1,-1;-2,1,1]*C_normed';
               alpha=std(S(1,:))/std(S(2,:));
               P=S(1,:)+alpha*S(2,:);
               l=size(block,1);
               H(1,f:f+l-1)=H(1,f:f+l-1)+(P-mean(P))./std(P);
            end
            pulse=double(H);
       end
   end
end

