% This code is developed by Adib Ahmadzadegan 
% all rights reserved.
% All rights reserved, please cite :
% 1)	Ahmadzadegan, A., Ardekani, A. M., & Vlachos, P. P. (2020). Estimation of the probability density function of random displacements from images. Physical Review E, 102(3), 033305.
% 2)	Ahmadzadegan, A., Ardekani, A., Vlachos, P. (2020). Estimation of the Probability Density Function of Random Displacements from Images. Purdue University Research Repository. doi:10.4231/34TJ-S109
% 3)	Ahmadzadegan, A., Harsa Mitra, Pavlos P. Vlachos, and Arezoo M. Ardekani. "Particle Image micro-Rheology (PIR) using displacement probability density function." Journal of Rheology 67, no. 4 (2023): 823-823.
% 4)	Ahmadzadegan, A., Sayantan Bhattacharya, Arezoo M. Ardekani, and Pavlos P. Vlachos. "Uncertainty estimation for ensemble particle image velocimetry." Measurement Science and Technology 33, no. 8 (2022): 085302.


clc
clear
close all


%% image direcctory information 
imdir='/scratch/shannon/a/aether/Projects/Microscale_measurments/micro-rheology/MonteCarlo/';%
N=2000 %number of frames
mintimelag=1
maxtimelag=25   % maximum correlation step 
stradlemod=0;  % 0 if the capture is continues and 1 if the images are captured through a stadling signal
%image file name base
filename=strcat('im_','Dt_',num2str(Dt),'_dp',num2str(dp(1)),'_');
numdig= '%01i'
filetype= '.png';

%% iPED settings :
dofilter=1;   % Edge tapering filter 
gausblur=0;   % gaussian blue for subresolved particles 
mbf=0;        % 1 for mean intensity subtraction 
maf=0;        % 1 for mean intensity subtraction 
ensemf=1;     % 1 for ensemble in Fourier  
windowz=64;   % correlation window size normally 64
mean=1;       % Mean image subtraction to remove the background 
pcolor=1;     % if 1 means particles are white, 0 means particles are black
outlier =0;   % 1 for outlier removal from images 
normalalize=0;% 1 for Intnsity normalization in images 
dosave=0;     % save preprocessed images. (Not required)

resave=strcat(resdir,'/results5000_',num2str(imstart));  % saveing directory
mkdir(resave); 
processed_dir=strcat(resdir,'/pre_processed/'); %directory for preprocessed images
mkdir(processed_dir)

%% estimate the background image
[bkga]=Meansub_v2(imdir,filename,numdig,filetype,processed_dir,mean,pcolor,outlier,N,...
    normalalize,dosave,imstart)
%% run iPED for different timelags
parfor corstep = mintimelag:maxtimelag
    MSD_cor(imdir,filename,numdig,filetype,resave,dofilter,...
        N,corstep,corstep,gausblur,mbf,maf,ensemf,windowz,stradlemod,imstart,bkga)
end


%%
addpath 'PIR\microrheology-master\src\Moduli';
addpath 'PIR\microrheology-master\src\MSD';

scale=3.5/20*10^-6; %pixel size in [m]
ms2=ms2*(scale)^2;%[m^2]
dt=1/10;% interframe time[s]
r=0.5*500*10^-9  % [m] Particle radius

timelags=(mintimelag:maxtimelag).*dt;       
msd=ms2(mintimelag:maxtimelag); 
Kb=1.38064852e-23;% m2 kg s-2 K-1
T=22.5 +273.15; %K   temperature of the system

removeoutlier=0;
matout=1;
%% estimate MSD from PDF
[msd]=PDF2MSD(removeoutlier,matout,resave,filename,mintimelag,maxtimelag,N);

%% making_logarithmically_spaced_msd_vs_tau
[msdtau] = making_logarithmically_spaced_msd_vs_tau(msd*10^12,timelags,max(timelags));
%% find complex viscosity Gp and Gpp as a function of frequency omega
[ omega,Gs,Gp,Gpp, dd, dda ] = calc_G(msdtau(:,1),msdtau(:,2),r*10^6,2,T,0.03,2);







