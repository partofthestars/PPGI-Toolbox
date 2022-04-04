# PPGI-Toolbox 
<p align="center"><img width=60% src="https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/Logo.jpg"></p>
<p align="center"><b>A MATLAB Toolbox for Photoplethysmography Imaging</b></p>
<p align="center">by Christian S. Pilz,<br>
  Aachen, 2019</p>
<p align="center">Version beta0.1</p>
<p align="center"><img width=10% src="https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/beta.jpg"></p>

## News

**04.04.2022**<br>
I uploaded the full implementation of the diffusion process model.<br>
An example implementation of a stochastic oscillator is added as well:
<p align="center"><img width=60% src="https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/results/oscillator.png"></p>

**29.11.2021**<br>
Our work was recently featured in the following Nature article about digital medicine (npj Digit. Med.):<br>
*Dasari, A., Prakash, S.K.A., Jeni, L.A. et al.. Evaluation of biases in remote photoplethysmography methods. npj Digit. Med. 4, 91 (2021).*  
[https://doi.org/10.1038/s41746-021-00462-z](https://doi.org/10.1038/s41746-021-00462-z)

## Supported by

| [![VideoBlocks](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/cancontrols.png)](http://cancontrols.com)  | [![GraphicStock](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/IXP.png)](http://ixp-duesseldorf.de) | [![AudioBlocks](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/UKAachen.png)](https://www.ukaachen.de/kliniken-institute/klinik-fuer-psychiatrie-psychotherapie-und-psychosomatik/forschung/neuropsychologie.html) | [![AudioBlockss](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/MedIT_RWTH.png)](https://www.medit.hia.rwth-aachen.de) | 
|:---:|:---:|:---:| :---:|

## Funded by

| [![BMBF](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/efre_eu_logo.png)](https://ec.europa.eu/regional_policy/en/funding/erdf/)  | [![EFRE](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/BMBF.png)](http://bmbf.de) | [![EFRE.NRW](https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/efre_nrw_logo.jpg)](http://efre.nrw.de) |
|:---:|:---:|:---:|


## Introduction

During the last years measuring blood volume changes and heart rate measurements from facial images
gained attention at top computer vision conferences frequently. Most of these contributions focus
on how to cope with motion like head pose variations and facial expressions since any kind of motion on a speciﬁc
skin region of interest will destroy the raw signal in a way that no reliable information can be extracted anymore.
Beside from being able to estimate vitality parameters like heart rate and respiration, the functional survey of wounds
as well as quantiﬁcation of allergic skin reaction are further topics of discovered employments of skin blood
perfusion analysis. Recently, prediction of emotional states, stress, fatigue and sickness became interesting new achievements in this area, pushing the focus of this technology further towards human-machine interaction.
In contrast to the genuine medical use-case of the technology, in computer vision and human-machine interaction we can’t expect any cooperative behavior of the user without introducing lack of convenience and a reduction of
the general user acceptance. Further, beyond any well tempered clinical and laboratory like scenarios, the majority
application will face strong challenging environmental changes and differences much more quite common. Thus,
there’s an emerging demand to produce better features and models signiﬁcant more robust to nuisance factors, still
preserving the desired target information. To reach such a formulation a fundamental profound understanding of the
underlying optical and mathematical properties is one of the current foci of this research discipline.
<br><br>
<p align="center"><img width=75% src="https://github.com/partofthestars/PPGI-Toolbox/blob/master/media/ico/cvpr2018.jpg"></p>
<i><b>Figure 1.</b> Rigid and non-rigid facial motions act as nuisance factors
on the tiny blood volume changes inherently destroying the target
information of heart rate under the conventional formulation of the
problem. Utilizing features invariant with respect to the action of the group of nuisance
transformations making it possible to estimate heart rate information under everyday facial motions. As illustrated for the above head motions, the green channel information doesn’t yield to reasonable heart rate information in the frequency domain. However, in case we're able to extract invariant features we obtain a clear signal in both time
and frequency domain.</i>



## Example Data

The example data can be download by the following link:<br>
https://www.dropbox.com/s/vv6ethy5az16wt4/example_data.mat?dl=0 <br>
(NOTE: Leave me a message in case the file isn't available anymore via dropbox!)<br>
Alternative download location: https://gw.cancontrols.com/LGI_DATABASE/toolbox/example_data.mat

Place the example_data.mat file into ./media/data/ folder.<br>
The example_data.mat contains the reference finger pulse oximeter waveform (PPG)
and the color image data of a face finder detection result in the rgb cell.

## Tests

In the tests directory there're separate example Matlab scripts. One for each algorithm respectively.
Initially the startup.m script must be executed to add all needed directories.
As second step executing test_skin.m will perform skin segmentation on the example_data.mat in case it was
downloaded and placed into the ./media/data/ folder. The obtained skin pixels will be stored into the example_data.mat file
and can be reused to test all other algorithm scripts.

## Algorithms

- Channel Mean (G) [8,9]
- Spatial Subspace Rotation (SSR) [6]
- Algorithmic Principles of Remote PPG (POS)[5]
- Local Group Invariance (LGI) [1,2]
- Diffusion Process (DP) [3,4]
- Riemannian-PPGI (SPH) [1]

## Evaluation

- Correlation
- Bland-Altman
- RMSE/ MSE
- SNR [7]

## Databases

- UBFC-RPPG [10]: [Download link](https://sites.google.com/view/ybenezeth/ubfcrppg)
- LGI Multi Session [1,2,3]: [Download link](https://github.com/partofthestars/LGI-PPGI-DB)

## References

1. Christian S. Pilz, Ibtissem Ben Makhlouf, Vladimir Blazek, Steffen Leonhardt, 
On the Vector Space in Photoplethysmography Imaging, 
The IEEE International Conference on Computer Vision (ICCV) Workshops,
Seoul/ South Korea, 2019
2. Christian S. Pilz, S. Zaunseder, J. Krajewski, V. Blazek,
Local Group Invariance for Heart Rate Estimation from Face Videos in the Wild,
The IEEE Conference on Computer Vision and Pattern Recognition (CVPR) Workshops, pp.1254-1262, Salt Lake City, 2018
3. Christian S. Pilz, Jarek Krajewski, Vladimir Blazek.
On the Diffusion Process for Heart Rate Estimation from Face Videos under Realistic Conditions.
Pattern Recognition: 39th German Conference, GCPR 2017, Basel, Switzerland.
Proceedings (Lecture Notes in Computer Science), pp. 361-373, Springer, 2017
4. Christian S. Pilz, Sebastian Zaunseder, Ulrich Canzler, Jarek Krajewski.
Heart rate from face videos under realistic conditions for advanced driver monitoring. 
Current Directions in Biomedical Engineering, De Gruyter, Berlin, pp. 483–487, 2017.
5. Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG. IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491
6. W. Wang, S. Stuijk and G. de Haan, "A Novel Algorithm for Remote Photoplethysmography: Spatial Subspace Rotation," in IEEE Transactions on Biomedical Engineering, vol. 63, no. 9, pp. 1974-1984, Sept. 2016.
7. De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886
8. Verkruysse, W., Svaasand, L. O., & Nelson, J. S. (2008). Remote plethysmographic imaging using ambient light. Optics express, 16(26), 21434-21445
9. M. Hülsbusch. A functional imaging technique for opto-electronic assessment of skin perfusion. PhD thesis, RWTH Aachen University, 2008.
10. S. Bobbia, R. Macwan, Y. Benezeth, A. Mansouri, J. Dubois, "Unsupervised skin tissue segmentation for remote photoplethysmography", Pattern Recognition Letters, 2017.

## License
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://opensource.org/licenses/gpl-3.0.html)

If you use this toolbox, please cite this paper:<br>
*Christian S. Pilz, Ibtissem Ben Makhlouf, Vladimir Blazek, Steffen Leonhardt, 
On the Vector Space in Photoplethysmography Imaging, 
The IEEE International Conference on Computer Vision (ICCV) Workshops,
Seoul/ South Korea, 2019*

## Results

Evaluation was conducted on the following databases:
- UBFC-RPPG
- LGI Multi Session

### Overview

Comparison of heart rate prediction accuracy utilizing different feature operators on diverse face video
databases. Each corresponding table entry represents the average Pearson’s correlation coefﬁcient together with the
average root-mean-square error (RMSE) value in BPM.

| Database  | Green [7,8] | SSR [5] | POS [3]| LGI [2] | SPH [1] |
| ------------- | ------------- |------------- |------------- |------------- |------------- |
| UBFC | 0.16/22.1 | 0.54/4.95 | 0.68/4.42 | 0.75/5.94 | 0.73/3.21 |
| LGI Resting | 0.41/2.61 | 0.49/1.99 | 0.41/2.10 | 0.69/1.41 | 0.71/1.49 | 
| LGI Rotation | 0.15/13.2 | 0.06/10.9 | 0.12/5.32 | 0.67/1.92 | 0.56/2.54 | 
| LGI Talk | 0.15/46.6 | 0.12/27.8 | 0.01/37.7 | 0.51/14.7 | 0.23/27.8 | 
| LGI Gym | 0.01/33.5 | 0.03/21.2 | 0.15/12.2 | 0.42/2.65 | 0.26/3.65 | 

In the following the box plot statistics for the UBFC and the LGI database are visualized. In addition to previous table the
inﬂuence of the Diffusion process [3,4] incorporated into the LGI and SPH approach is constituted. Both approaches beneﬁt
from its stochastic interpretation as quasi-periodic nature of blood volume changes in contrast to the traditional standard Fourier based spectral peak search.

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
