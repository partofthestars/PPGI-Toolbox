# PPGI-Toolbox 
<p align="center"><img width=60% src="https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/Logo.jpg"></p>
<p align="center"><b>A MATLAB Toolbox for Photoplethysmography Imaging</b></p>
<p align="center">by Christian S. Pilz,<br>
  Aachen, 2019</p>
<p align="center">Version beta0.1</p>
<p align="center"><img width=10% src="https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/beta.jpg"></p>

## Supported by

| [![VideoBlocks](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/cancontrols.png)](http://cancontrols.com)  | [![AudioBlocks](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/MedIT_RWTH.png)](https://www.medit.hia.rwth-aachen.de) | [![GraphicStock](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/IXP.png)](http://ixp-duesseldorf.de) |
|:---:|:---:|:---:|


## Funded by

| [![BMBF](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/efre_eu_logo.png)](https://ec.europa.eu/regional_policy/en/funding/erdf/)  | [![EFRE](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/BMBF.png)](http://bmbf.de) | [![EFRE.NRW](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/efre_nrw_logo.jpg)](http://efre.nrw.de) |
|:---:|:---:|:---:|


## Example Data

The example data can be download by the following link:<br>
https://www.dropbox.com/s/vv6ethy5az16wt4/example_data.mat?dl=0

Place the example_data.mat file into ./media/data/ folder.<br>
The example_data.mat contains the reference finger pulse oximeter waveform (ppg)
and the color image data of a face finder detection result in the rgb cell.

## Algorithms

- Channel Mean (G) [7,8]
- Spatial Subspace Rotation (SSR) [5]
- Algorithmic Principles of Remote PPG (POS)[3]
- Local Group Invariance (LGI) [1,2]
- Diffusion Process (DP) [4]
- Riemannian-PPGI (SH) [1]

## Evaluation

- Correlation
- Bland-Altman
- RMSE/ MSE
- SNR [6]

## Databases

- UBFC-RPPG [9]: [Download link](https://sites.google.com/view/ybenezeth/ubfcrppg)
- LGI Multi Session [1,2,4]: [Download link](https://github.com/partofthestars/LGI-PPGI-DB)

## References

1. Christian S. Pilz, Vladimir Blazek, Steffen Leonhardt, On the Vector Space in Photoplethysmography Imaging,<br>
Preprint: arXiv:1903.03316 [cs.CV], 2019
2. Christian S. Pilz, S. Zaunseder, J. Krajewski, V. Blazek,
Local Group Invariance for Heart Rate Estimation from Face Videos in the Wild,
The IEEE Conference on Computer Vision and Pattern Recognition (CVPR) Workshops, pp.1254-1262, Salt Lake City, 2018
3. Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG. IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491
4. Christian S. Pilz, Jarek Krajewski, Vladimir Blazek.
On the Diffusion Process for Heart Rate Estimation from Face Videos under Realistic Conditions.
Pattern Recognition: 39th German Conference, GCPR 2017, Basel, Switzerland.
Proceedings (Lecture Notes in Computer Science), Springer, 2017
5. W. Wang, S. Stuijk and G. de Haan, "A Novel Algorithm for Remote Photoplethysmography: Spatial Subspace Rotation," in IEEE Transactions on Biomedical Engineering, vol. 63, no. 9, pp. 1974-1984, Sept. 2016.
6. De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886
7. Verkruysse, W., Svaasand, L. O., & Nelson, J. S. (2008). Remote plethysmographic imaging using ambient light. Optics express, 16(26), 21434-21445
8. M. Hülsbusch. A functional imaging technique for opto-electronic assessment of skin perfusion. PhD thesis, RWTH Aachen University, 2008.
9. S. Bobbia, R. Macwan, Y. Benezeth, A. Mansouri, J. Dubois, "Unsupervised skin tissue segmentation for remote photoplethysmography", Pattern Recognition Letters, 2017.

## License
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://opensource.org/licenses/gpl-3.0.html)

If you use this dataset, please cite this paper:<br>
Christian S. Pilz, Vladimir Blazek, Steffen Leonhardt, On the Vector Space in Photoplethysmography Imaging, arXiv:1903.03316 [cs.CV], 2019

## Results

Vealuation was conducted on the following databses:
- UBFC-RPPG
- LGI Multi Session

### Overview

Comparison of heart rate prediction accuracy utilizing different feature operators on diverse face video
databases. Each corresponding table entry represents the average Pearson’s correlation coefﬁcient together with the
average root-mean-square error (RMSE) value in BPM.

| Database  | Green [7,8] | SSR [5] | POS [3]| LGI [2] | SPH [1] |
| ------------- | ------------- |------------- |------------- |------------- |------------- |
| Content Cell  | Content Cell  | Content Cell  | Content Cell  | Content Cell  | Content Cell  |
| Content Cell  | Content Cell  | Content Cell  | Content Cell  | Content Cell  | Content Cell  |


In the following the box plot statistics for the UBFC and the LGI database are visualized. In addition to previous table the
inﬂuence of the Diffusion process [4] incorporated into the LGI and SPH approach is constituted. Both approaches beneﬁt
from its stochastic interpretation as quasi-periodic nature of blood volume changes.

#### UBFC-RPPG:

| Pearson             |  RMSE |
:-------------------------:|:-------------------------:
![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/UBFC/ubfc_pearson.png)  |  ![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/UBFC/ubfc_rmse.png)

#### LGI Multi Session:

<br>
- Session 1: Head Resting

| Pearson             |  RMSE |
:-------------------------:|:-------------------------:
![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_Office_Head_Resting_pearson.png)  |  ![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_Office_Head_Resting_rmse.png)

<br>
- Session 2: Head Rotation


| Pearson             |  RMSE |
:-------------------------:|:-------------------------:
![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_Head_Rotation_pearson.png)  |  ![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_Head_Rotation_rmse.png)

<br>
- Session 3: Bicycle Ergometer


| Pearson             |  RMSE |
:-------------------------:|:-------------------------:
![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_GYM_pearson.png)  |  ![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_GYM_rmse.png)

<br>
- Session 4: Outdoor City Talk
<br>

| Pearson             |  RMSE |
:-------------------------:|:-------------------------:
![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_City_pearson.png)  |  ![](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/LGI/LGI_City_rmse.png)
