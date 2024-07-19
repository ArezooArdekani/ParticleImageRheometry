This document and the codes in this folder are written by Adib Ahmadzadegan, aahmadza@purdue.edu
All rights reserved, please cite :
1)	Ahmadzadegan, A., Ardekani, A. M., & Vlachos, P. P. (2020). Estimation of the probability density function of random displacements from images. Physical Review E, 102(3), 033305.
2)	Ahmadzadegan, A., Ardekani, A., Vlachos, P. (2020). Estimation of the Probability Density Function of Random Displacements from Images. Purdue University Research Repository. doi:10.4231/34TJ-S109
3)	Ahmadzadegan, Adib, Harsa Mitra, Pavlos P. Vlachos, and Arezoo M. Ardekani. "Particle Image micro-Rheology (PIR) using displacement probability density function." Journal of Rheology 67, no. 4 (2023): 823-823.
4)	Ahmadzadegan, Adib, Sayantan Bhattacharya, Arezoo M. Ardekani, and Pavlos P. Vlachos. "Uncertainty estimation for ensemble particle image velocimetry." Measurement Science and Technology 33, no. 8 (2022): 085302.
To run the code, open Matlab and go to the directory where the code is stored. 
Open the Demo code, that covers the steps to the full microrheology analysis using PIR. 
The main functions are Meansub_v2,MSD_Cor, PDF2MSD,and the microrheology master package.

The Demo code will first find the background image and subtract that from the images in the MSD_cor function. MSD cor is basically the iPED calculation for each time lag. Once all pdfs are calculated and saved, PDF2MSD uses those to calculate the MSD. Then MSD is used to find G’ and G” using the microrheology master package. All of these steps are included in the Demo code.
The outputs will be saved in your image folder. 
