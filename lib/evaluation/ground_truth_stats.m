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
%	ground_truth_stats.m:
%
% Description:
%
%   A static MATLAB class to compute ground truth statistics
%

classdef ground_truth_stats
    
   methods(Static)
       function [pearson, rmse, snr, snr_std, bpm] = get(ppg,waveform, fs)
           y=ppg;
           L=length(y);
           NFFT = 2^nextpow2(L);
           f = fs/2*linspace(0,1,NFFT/2+1);

           [spec_ppg]=spectrogram(y,256,240,f,fs,'yaxis');

           [peak_values, peak_frequencies] = max(abs(spec_ppg));
           bpm_ppg=f(peak_frequencies)*60;

           y=waveform;
           L=length(y);
           NFFT = 2^nextpow2(L);
           f = fs/2*linspace(0,1,NFFT/2+1);

           [spec,freqs,tmp, pxx]=spectrogram(y,256,240,f,fs,'yaxis');
           
         
           [peak_values, peak_frequencies] = max(abs(spec));
           bpm=f(peak_frequencies)*60;

           num_elements=size(bpm_ppg,2);
           if num_elements>size(bpm,2)
               num_elements=size(bpm,2);
           end

           [pearson]=corr(bpm_ppg(1,1:num_elements)',bpm(1,1:num_elements)');
           [rmse]=sum(sqrt((bpm_ppg(1,1:num_elements)-bpm(1,1:num_elements)).^2))/num_elements;
       
           %snr
           hr_f=bpm_ppg/60;
           for i=1:num_elements
               GTMask1 = (freqs >= hr_f(1,i)-0.1)&(freqs <= hr_f(1,i)+0.1);
               GTMask2 = (freqs >= hr_f(1,i)*2-0.2)&(freqs <= hr_f(1,i)*2+0.2);
               %SPower = sum(abs(spec(GTMask1|GTMask2,i)));
              
               %SPower = sum(abs(spec(GTMask1,i)));
               signal_power = sum(pxx(GTMask1,i));
               FMask2 = (freqs >= 0.5)&(freqs <= 4);
               %total_power = sum(abs(spec(FMask2,i)));
               total_power = sum(pxx(FMask2,i));
               SNR(i) = pow2db(signal_power/(total_power-signal_power));
           end
          
           [snr] = mean(SNR);
           [snr_std] =std(SNR);
       end
   end
end

