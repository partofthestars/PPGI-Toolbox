
# PPGI-Toolbox
<b>A MATLAB toolbox for Photoplethysmography Imaging</b><br>
<br>
by Christian S. Pilz, Aachen, 2019
<br>

## Example Data

The example data can be download given the following link:<br>
https://www.dropbox.com/s/vv6ethy5az16wt4/example_data.mat?dl=0

Place the example_data.mat file into this folder.
example_data.mat contains the reference finger pulse oximeter waveform (ppg)
and the rgb image data of a face finder detection result.

## Test skin

This will perfom simple skin segmentation on the sample rgb data of the example_data.mat.
The script will store the segmented skin pixels into the example_data.mat.
The skin pixels cell array can be used for the other test scripts without the need for
recomputing the skin segmentation again.
