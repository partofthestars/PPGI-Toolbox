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
%	separate.m:
%
% Description:
%
%   Implements the frequency estimation based upon the diffusion process algorithm.
%
% References:
%
%	Christian S. Pilz, Jarek Krajewski, Vladimir Blazek.
%   On the Diffusion Process for Heart Rate Estimation from Face Videos under Realistic Conditions.
%   Pattern Recognition: 39th German Conference, GCPR 2017, Basel, Switzerland.
%   Proceedings (Lecture Notes in Computer Science), Springer, 2017
%

classdef diffusion_process
    
   methods
       
       function obj = diffusion_process()
       
       end
       
       function [trace, obj] = get(obj, Y,dt,freqlist,nharm,BQ,sd,qr,ptrans,poverall)
           
            %sd = 0.025;
            %BQ = 0.001;   

            %qr = 0.0001;
            %Nimm=1;
            %dt=1/fs;
            %ptrans = 0.15;
            %poverall = 0;
            %freqlist=45:120;

            BF = [0 1; 0 0];  
            BL = [0;1]; 
            BH = [1 0];

            [frequencyResults] = track( ...
                                         Y,   ...
                                         dt,          ...
                                         freqlist,    ...
                                         nharm,        ...
                                         BF,          ...
                                         BQ,          ...
                                         BL,          ...
                                         BH,          ...
                                         sd^2,        ...
                                         qr,          ...
                                         ptrans,      ...
                                         poverall);
             
            [S,CS,QCS] = separate( ...
                                     Y, ...
                                     dt,      ...
                                     frequencyResults,    ...
                                     nharm,        ...
                                     BF,          ...
                                     BQ,          ...
                                     BL,          ...
                                     BH,          ...
                                     sd^2,        ...
                                     qr);


             trace=CS;
       end
   end 
end

