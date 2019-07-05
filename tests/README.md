
# PPGI-Toolbox
<b>A MATLAB toolbox for Photoplethysmography Imaging</b><br>
<br>
by Christian S. Pilz, Aachen, 2019<br>
<br>

## Example Data

The example data can be download given the following link:<br>
https://www.dropbox.com/s/vv6ethy5az16wt4/example_data.mat?dl=0

Place the example_data.mat file into this folder.
example_data.mat contains the reference finger pulse oximeter waveform (ppg)
and the rgb image data of a face finder detection result.

## Setup
Execute startup.m from the toolbox root dir to set all path needed by the scripts

## Test skin

This will perfom simple skin segmentation on the sample rgb data of the example_data.mat.
The script will store the segmented skin pixels into the example_data.mat.
The skin pixels cell array can be used for the other test scripts in order 
to avoid the need for a recomputation of the skin segmentation.
this is essentially pretty usefull during extensive data evaluation procedures
where often several attemps using different parameters are needed to obtained
a good calibrated result.

## Test diffusion process

This will track and separate the pulse signal. The important parameters of the model are all explained in the script.
Dependent on the contamination of the source signal the parameters must be adjusted to obtain good results.
Unfortunately, i wasn't able to find a closed form solution in determining proper parameter settings. I spent a lot of time
to collect the necessary experience about their influence and how to set these properly.
