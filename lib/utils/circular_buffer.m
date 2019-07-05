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
%	circular_buffer.m:
%
% Description:
%
%   A  MATLAB circular buffer class to handle streaming data
%

classdef circular_buffer
    
   properties
        buf_;
        head_;
        tail_;
        max_size_;
        full_;
        cell_
   end
   
   methods
       
       function obj = circular_buffer(max_size,cell)
           
          if nargin > 0
             if isnumeric(max_size)
                obj.max_size_ = max_size;
                obj.buf_= [];
                obj.head_= 1;
                obj.tail_= 1;
                obj.full_= logical(false);
                obj.cell_= logical(cell);
             else
                error('size must be numeric')
             end
          end
          
       end
       
       function ret = capacity(obj)
           ret = obj.max_size_;
       end
       
       function ret =size(obj)
           ret = obj.max_size_;

            if obj.full_== 0
                if obj.head_ >= obj.tail_
                    ret = obj.head_ - obj.tail_;
                else
                    %ret = obj.max_size_ + obj.head_ - obj.tail_;
                    ret = obj.head_ - obj.tail_;
                end
            end
       end
       
       function obj = put(obj,data)
           if obj.cell_
               obj.buf_{obj.head_} = data;
           else
               obj.buf_(obj.head_,:) = data;
           end
          
%            if obj.full_
%                obj.tail_ = mod(obj.head_, obj.max_size_)+1; 
%            end
           
           obj.head_ = mod(obj.head_, obj.max_size_)+1;
           obj.full_ = (obj.head_ == obj.tail_);
       end
       
       function [values, obj] = get(obj,window_size,shift)
           if obj.empty()
               %values=[];
               return;
           end
           try
               if obj.cell_
                   values=obj.buf_(1,obj.tail_:obj.tail_+window_size-1)
               else
                   if (obj.tail_+window_size-1)<=obj.max_size_
                       values=obj.buf_(obj.tail_:obj.tail_+window_size-1,:);
                   else
                       tail_delta=obj.max_size_-obj.tail_;
                       values=obj.buf_(obj.tail_:obj.tail_+tail_delta,:);
                       head_delta=window_size-tail_delta-1;
                       values=[values; obj.buf_(1:head_delta,:)];
                   end
               end
           catch exception
                throw(exception)       
           end
           
           obj.tail_ = mod(obj.tail_ + shift, obj.max_size_);
           if obj.tail_ == 0
               obj.tail_=1;
           end
       end
       
       function obj = reset(obj)
           obj.head_ = obj.tail_;
           obj.full_ = logical(false);
       end
  
       function ret = empty(obj)
           %if head and tail are equal, we are empty
           ret = (obj.full_ == 0) & (obj.head_ == obj.tail_);
       end
       
       function ret = full(obj)
           %if tail is ahead the head by 1, we are full
           ret = obj.full_;
       end
   end
end

