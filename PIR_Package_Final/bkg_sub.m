%Background subtraction mean/min sub method
%This code imports images and find the mean/min image and then substracts that
%from every image.
%This code was written by Adib Ahmadzadegan 1/24/18
%

clear
clc
%imdir='E:/PID/phantom/590nm/pre_processed/';  % on windows use /
imdir='E:/bacteria_data/ecoli/pdf_test/Inverted_bkgsub/';
im_file_name='pdf_test_%04d.tif';
%im_file_name='even_%04d.tif'; % image file name
processed_dir=strcat(imdir,'preprocessed/')
mkdir(processed_dir)

mean=1;%%if you want mean subtraction set mean to 1
pcolor=1;%%what color is your signal ? it the particles are white set pcolor=1 else
%%if particle is black pcolor=0
N_image=999; %number of images to be processed

outlier =0 ; % if in your images there is black or white broken pixel which you
%%want to change to background color
normalalize=0; % set normalize=1 if you want to normalize the image afterward.
%%It will brighten the signal


flag =0;
if mean ==1
    flag = 1;
end
sprintf(strcat(imdir,im_file_name),1)
image =imread(sprintf(strcat(imdir,im_file_name),1));
imagesize_x=size(image,1);
imagesize_y=size(image,2);


if outlier==1
    readim=sprintf(strcat(processed_dir,'even_bkg/'));
    mkdir(readim)
end


h = waitbar(0,'calculating...','Name','Please wait...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for kstep=0:100:0
    bkg=zeros(imagesize_x,imagesize_y);
    minim=2^16*ones(imagesize_x,imagesize_y);
    
    
    for jj=1:N_image
        %         h =waitbar(jj/N_image,h,strcat('Mean Image', {' '},num2str(jj),' out of',{' '},num2str(N_image), ' is done. '))
        image =imread(sprintf(strcat(imdir,im_file_name),jj+kstep));
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
            mm=mean2(IM2);
            for x=1:imagesize_x
                for y=1:imagesize_y
                    if IM2(x,y)<mm
                        IM2(x,y)=0;%1.5*mu;
                    end
                end
            end
            fname=sprintf(strcat(readim,'even_bkg_','%d.tif'),jj+kstep);
            imwrite(uint16(IM2),fname);
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
    
    delete(h)
    if flag==1
        bkga=bkg./N_image;  % avereged bkg
        bkgdir=strcat(processed_dir,'bkg',num2str(kstep),'.tif');
        imwrite(uint16(bkga),bkgdir);
    else
        bkgdir=strcat(processed_dir,'min',num2str(kstep),'.tif');
        imwrite(uint16(minim),bkgdir);
    end
    h = waitbar(0,'calculating...','Name','Please wait...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    for t=1:N_image
        h =waitbar(t/N_image,h,strcat('Image MeanSub', {' '},num2str(t),' out of',{' '},num2str(N_image), ' is done. '));
        
        if outlier==1
            fname=sprintf(strcat(readim,'even_bkg_','%d.tif'),t+kstep);
            image =(imread(fname));
        else
            image =imread(sprintf(strcat(imdir,im_file_name),t+kstep));
        end
        
        if pcolor==1
            IM2 =image;
        else
            IM2 =imcomplement(image);
        end
        
        IM2 = double(((IM2)));
        if flag==1
            fimage = IM2-(bkga);
            fname=sprintf(strcat(processed_dir,'mean_enh_','%d.tif'),t+kstep)
        else
            fimage = IM2-(minim);
            fname=sprintf(strcat(processed_dir,'min_enh_','%d.tif'),t+kstep)
        end
        
        if normalalize==1
            fimage_f=(fimage./max(max(fimage))).*2^16;
        else
            fimage_f=fimage;
        end
        imwrite(uint16(fimage_f),fname);
    end
    
    delete(h)
    
end




