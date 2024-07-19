function [bkga]=Meansub_v2(imdir,filename,numdig,filetype,processed_dir,mean,pcolor,outlier,N_image,...
    normalalize,dosave,imstart)

% %Background subtraction mean/min sub method
% %This code imports images and find the mean/min image and then substracts that
% %from every image.
% %This code was written by Adib Ahmadzadegan 1/24/18
% %
%
% imdir='C:/Users/Adib/Desktop/Projects/PIR/1micron/raw/';
% im_file_name=strcat(filename,numdig,filetype)%'1micron_water_diffusion_60X_40pps_40microsecexpo_%04d.tif'; % image file name
% (strcat(im_file_name,1))
% strcat(sprintf(imdir,im_file_name),1)



% processed_dir=strcat(imdir,'pre_processed/')
% mkdir(processed_dir)
%
%
% mean=0;%%if you want mean subtraction set mean to 1
% pcolor=1;%%what color is your signal ? it the particles are white set pcolor=1 else
% %%if particle is black pcolor=0
% N_image=100; %number of images to be processed
% outlier =1 ; % if in your images there is black or white broken pixel which you
% %%want to change to background color
% normalalize=0; % set normalize=1 if you want to normalize the image afterward.
% %%It will brighten the signal


flag =0;
if mean ==1
    flag = 1;
end



filenum=num2str(1,numdig);
matFilename  = strcat(filename,filenum,filetype);
image =imread([imdir, matFilename]);%imread(strcat(sprintf(imdir,im_file_name),1));

imagesize_x=size(image,1);
imagesize_y=size(image,2);
bkg=zeros(imagesize_x,imagesize_y);
minim=(2^16)*ones(imagesize_x,imagesize_y);

if outlier==1
    readim=sprintf(strcat(processed_dir,'even_bkg/'));
    mkdir(readim);
end



for jj=imstart:N_image+imstart-1
    filenum=num2str(jj,numdig)
    matFilename  = strcat(filename,filenum,filetype);
    image =imread([imdir, matFilename]);%imread(sprintf(strcat(imdir,im_file_name),j));
            if size(image,3)==3
%             disp('rgb')
            image=rgb2gray(image);
        end
    if pcolor==1
        IM2 =image;
    else
        IM2 =imcomplement(image);
    end
    
    
    if outlier==1
        % Compute the Otsu level
        [counts,ss] = imhist(IM2);
        f = fit(ss,counts,'gauss1');
        mu=f.b1;
        sigma=f.c1/sqrt(2);
        IM2 = (double(IM2));
        for x=1:imagesize_x
            for y=1:imagesize_y
                if IM2(x,y)<mu
                    IM2(x,y)=mu;
                end
            end
        end
        
        
        filenum=num2str(jj,numdig);
        matFilename  = strcat('even_bkg_',filenum,filetype);
        fname =[readim, matFilename];
        imwrite(uint8(IM2),fname);
    end
    
    if flag==1
        bkg=bkg+double(IM2);
    else
        for i=1:size(IM2,1)
            for j=1:size(IM2,2)
                if minim(i,j)>IM2(i,j)
                    minim(i,j)=IM2(i,j);
                end
            end
        end
    end
end

if flag==1
    bkga=bkg./N_image;  % avereged bkg
    bkga2save=uint16(bkga);
    matFilename  = strcat('mean_image',filetype);
    fname =[processed_dir, matFilename];
    bkgadir=strcat(processed_dir,'bkga');
    save(bkgadir,'bkga')
    imwrite(bkga2save,fname);
else
    bkga=minim;
    matFilename  = strcat('min_image',filetype);
    fname =[processed_dir, matFilename];
    imwrite(uint16(minim),fname);
    
end

if dosave==1
    for t=imstart:N_image+imstart-1
        if outlier==1
            filenum=num2str(t,numdig);
            matFilename  = strcat('even_bkg_',filenum,filetype);
            image =imread([imdir, matFilename]);
                       if size(image,3)==3
%             disp('rgb')
            image=rgb2gray(image);
        end
        else
            
            
            filenum=num2str(t,numdig);
            matFilename  = strcat(filename,filenum,filetype);
            image =imread([imdir, matFilename]);
                       if size(image,3)==3
%             disp('rgb')
            image=rgb2gray(image);
        end
            
        end
        
        if pcolor==1
            IM2 =image;
        else
            IM2 =imcomplement(image);
        end
        
        IM2 = double(((IM2)));
        if flag==1
            fimage = IM2-(bkga);
            
            filenum=num2str(t,numdig)
            matFilename  = strcat('mean_enh_',filenum,filetype);
            fname =[processed_dir, matFilename];
        else
            fimage = IM2-(minim);
            
            filenum=num2str(t,numdig)
            matFilename  = strcat('min_enh_',filenum,filetype);
            fname =[processed_dir, matFilename];
        end
        
        if normalalize==1
            fimage_f=(fimage./max(max(fimage))).*256;
        else
            fimage_f=fimage;
        end
        imwrite(uint16(fimage_f),fname);
    end
end
