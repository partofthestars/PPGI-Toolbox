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
%   A static MATLAB class to extract skin pixels based upon
%   YCbCr color space thresholding
%

classdef skin
    
   methods(Static)
       function out = get(rgb)
           %convert the image from RGB to YCbCr
           img_ycbcr = rgb2ycbcr(rgb);
           Cb = img_ycbcr(:,:,2);
           Cr = img_ycbcr(:,:,3);
           
           %detect skin
           [r,c,v] = find(Cb>=98 & Cb<=142 & Cr>=133 & Cr<=177);
            numind = size(r,1);
            
            %extract skin Pixels
            for i=1:numind
                out(i,:)=double(rgb(r(i),c(i),:));
            end
       end
   end
end

