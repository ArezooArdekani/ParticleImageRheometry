clear
clc
imdir='E:\PID\phantom\194nm\';
im_file_name='194nm_60x_50uinside_21c_';
numdig= '%04i';
filetype='.tif'; % image file name
processed_dir=strcat(imdir,'pre_processed/')
mkdir(processed_dir)
readim=sprintf(strcat(processed_dir,'even_bkg/'));
mkdir(readim);


for i=1:5000 
filenum=num2str(i,numdig);
matFilename  = strcat(im_file_name,filenum,filetype);
image =imread([imdir, matFilename]);%imread(strcat(sprintf(imdir,im_file_name),1));
image(image < mean2(image))= mean2(image);
matFilename  = strcat('even_',filenum,filetype);
fname =[processed_dir, matFilename];
imwrite(uint16(image),fname);
end