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
%	spatial_subspace_rotation.m:
%
% Description:
%
%   Implements the spatial subspace rotation algorithm from TU Eindhoven
%
% References:
%
%   W. Wang, S. Stuijk and G. de Haan, "A Novel Algorithm for Remote Photoplethysmography: Spatial Subspace Rotation," 
%   in IEEE Transactions on Biomedical Engineering, vol. 63, no. 9, pp. 1974-1984, Sept. 2016.
%

classdef spatial_subspace_rotation
    
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
       
       function obj = spatial_subspace_rotation(fs,window_size,overlap)
           
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
     
       function [pulse, obj] = get(obj,varargin)%skin_pixels)
           
           pulse=0;%[0 0 0];
           
           skin_pixels=varargin{1};

           [rows cols]=size(skin_pixels);

           %spatial RGB correlation:
           C=(skin_pixels'*skin_pixels)/(rows*cols);
           [V,D] = eig(C);

           [D,I] = sort(diag(D),'descend');
           V = V(:, I);

           obj.U_=obj.U_.put(V);
           obj.S_=obj.S_.put(D);

           if obj.U_.size()>1

               %rotation between the skin vector and orthonormal plane
               [U obj.U_]=obj.U_.get(2,1);
               R=[U{2}(:,1)'*U{1}(:,2) U{2}(:,1)'*U{1}(:,3)];
               %scale change
               [S obj.S_]=obj.S_.get(2,1);
               SR=[sqrt(S{2}(1)/S{1}(2)) sqrt(S{2}(1)/S{1}(3))].*R;
               obj.SR_B_=obj.SR_B_.put(SR*[U{1}(:,2)'; U{1}(:,3)']); 

               if obj.SR_B_.size()>obj.fs_*(obj.window_size_)+obj.shift_
                   [frame obj.SR_B_]=obj.SR_B_.get(obj.fs_*obj.window_size_,obj.shift_);
                   sigma=std(frame(:,1))/std(frame(:,2));  
                   tmp=frame(:,1)-sigma*frame(:,2);
                   pulse=tmp-mean(tmp);
                   pulse=pulse(end,1);
               end
           end              
       end
       
       %function obj = reset(obj)
       %    obj.U_.reset();
       %end
   end
end

