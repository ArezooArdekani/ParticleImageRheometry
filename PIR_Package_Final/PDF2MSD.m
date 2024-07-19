function [msd]=PDF2MSD(removeoutlier,matout,resave,filename,mint,maxt,N)
%removeoutlier=0; if you want to remove outliers in the PDF set this equals
%to 1
%matout=1;   If you want to save the results
%filedir=E:\imagefilesMP1_Dt0.5_runnumber_1\dofilter1_gaussb0_mbf1maf0ensemf1\;
% filename='pdfofensemble_Adib_999_sat100window64_'%'PDF_ensemble_5000_'%
% mint = minimum correlation step size
% maxt = maximum correlation step size
% N number of correlation ensembles

for i=mint:maxt
    filedir=strcat(resave,'\','pdf_iPED_',num2str(N-i),'window',num2str(64),'_','timelag',num2str(i)),'.mat');
    load(filedir);
    imsize=size(pdf_iPED,1);
    cropsize=imsize;%200
    [X,Y]=meshgrid(1:cropsize);
    
    %thresholding the PDF to clean up
    PDF=(pdf_iPED);
    thvalue=abs(min(PDF(:)));%abs(quantile(PDF(:),0.01));
    PDF(PDF<thvalue)=thvalue;
    PDF=PDF-thvalue;
    F=PDF;
    [F]= FITG(F);    % Dedrift the PDF by shifting it to the center using a 2D Gaussian fit
    figure(9)
    subplot(1,4,1)
    imagesc(pdf_iPED)
    subplot(1,4,2)
    mesh(pdf_iPED)
    figure(9)
    subplot(1,4,3)
    mesh(F)
    subplot(1,4,4)
    imagesc(F)
    F=F./trapz(F(:));    % make sure the PDF has a integral of 1 
    % estimate the MSD from the PDF: 
    
    MS=abs((X-cropsize/2-1).*(Y-cropsize/2-1).*F)*4;    
    msd(i)=trapz(MS(:));

end

savemsd  = strcat(resave,'\msd_N',num2str(N));
save(savemsd,'msd');

end

function [J,xcenter,ycenter]= FITG(Z)
Z=(Z);%-mean(Z(:));
[X,Y]=meshgrid(1:size(Z,2),1:size(Z,1));
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
MdataSize = min(size(X,1),size(Y,2));
lb = [0,0,0,0,0,-inf];
ub = [realmax('double'),MdataSize,(MdataSize)^2,MdataSize,(MdataSize)^2,inf];
x0 = [max(Z(:)),size(Z,2)/2,5,size(Z,1)/2,5,0];
opts = optimset('Display','off');
[x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub,opts);
XC=size(Z,2)/2+1;
YC=size(Z,2)/2+1;
ycenter=YC-x(4);
xcenter=XC-x(2);

J = imtranslate(Z,[xcenter,ycenter]);

end

