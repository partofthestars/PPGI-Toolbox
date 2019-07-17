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
%	spherical_mean.m:
%
% Description:
%
%   The spherical operator on the Riemann manifold: Riemannian PPGI.
%   NOTE: The current implementation is a rough approximation
%   of the papers mathematical formulation.
%
% References:
%
%   Christian S. Pilz, Vladimir Blazek, Steffen Leonhardt.
%   On the Vector Space in Photoplethysmography Imaging, 
%   Preprint: arXiv:1903.03316 [cs.CV], 2019
%

classdef spherical_mean
    
   methods
       
       function obj = spherical_mean()
       
       end
       
       function [mu, obj] = get(obj,rgb)
           
           sphere(1:size(rgb,1),1:3)=0;
           for i=1:size(rgb,1)
               sphere(i,:)=(rgb(i,:)./norm(rgb(i,:)));
           end
           
           expectation=mean(sphere);
           
           [mu(1,1),mu(1,2),r] = cart2sph(expectation(:,1),expectation(:,2),expectation(:,3));
       end
   end 
end

