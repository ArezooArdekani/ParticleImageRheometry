function []=MSD_cor(imdir,filename,numdig,filetype,resave,dofilter,N,tmin,tmax,gausblur,mbf,maf,ensemf,windowz,stradlemod,imstart,bkga)
%         imdir ='E:\bacteria_data\ecoli\pdf_test\';
%         filename='20x_hcb1736_secondaryculture_od0_t';
%         numdig= '%04i'
%         filetype= '.tif'
%         resave=strcat(imdir,'\results');
%         mkdir(resave);
%         N=100; numebr of images
%         tmax, tmin are min and max time lags
%         dofilter=1;    gaussian Filter on the whole image
%         gausblur=0; % Gaussian blur for unresolved particles
%         mbf=1; mean before Gaussian filter
%         maf=1; mean after Gaussianfilter
%         ensemf=1;  ensemble in Fourier
%         bkga 
tic
% 
% h = waitbar(0,'calculating...','Name','Please wait...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

% imstart=1;       % the number of the first image
Montecarlo=0;    % to save all the correlations
cropp=1;         % window cropping for
% windowz=1024;      % Window size
fittt=1;       % if you want to fit Gaussian function to the PDF
saturation=0;  % This is to saturate images. Don't use this, it was for developmental studies.
%satratio=0.1:0.1:1
si=10;
if saturation==0
    si=100;
end
l=0;filter=ones;
breakflag=0;
kahler=zeros;
pid=zeros;
straddling=stradlemod;  % 1 if frame1 can only be correlated to frame2
convergencecriteria=0;
converged=0;
FFFF=zeros;
eeerrr=zeros;

for i=tmin:tmax     %%for MSD correlations
    % define a whole bunch of parameters, names are self explanatory
    max_intensity=0;
    SCC_ensembler=zeros;
    AC_ensembler=zeros;
    SCC_ensemblerr=zeros;
    AC_ensemblerr=zeros;
    PDF_GCCr=zeros;
    PDF_Adibr=zeros;
    PDF_GCCrr=zeros;
    PDF_Adibrr=zeros;
    SCC_ensembleFr=zeros;
    AC_ensembleFr=zeros;
    SCC_ensembleFrr=zeros;
    AC_ensembleFrr=zeros;
    AC=zeros;
    gcc_counter=0;
    adib_counter=0;
    AC_ensemble_kahler=zeros;
    SCC_ensemble_kahler=zeros;
    meannoiselevel=zeros;
    integrallevel=zeros;
    L_infinity=zeros;    
    ensgrids=0;
    CorCounter=0;
    for j=imstart:1+straddling:imstart+N-i
        CorCounter=CorCounter+1;
        %Read image 1 anad 2
        filenum=num2str(j,numdig)
        matFilename  = strcat(filename,filenum,filetype);
        [xr1_0,map]=imread([imdir, matFilename]);
        if size(xr1_0,3)==3
%             disp('rgb')
            xr1_0=rgb2gray(xr1_0);
        end
        xr1_0_m=double(xr1_0)-double(bkga);
        xr1_0=double(xr1_0);
        filenum=num2str(j+i,numdig);
        matFilename  = strcat(filename,filenum,filetype);
        [xr2_0,map]=imread([imdir, matFilename]);
         if size(xr2_0,3)==3
            xr2_0=rgb2gray(xr2_0);
        end
        xr2_0_m=double(xr2_0)-double(bkga);
        xr2_0=double(xr2_0);
        
        kmax=size(xr1_0_m);

        for k=1:windowz:kmax(2)-windowz+1
            k;
            for kk=1:windowz:kmax(1)-windowz+1
                %                 if cropp==1
                xr1_m=imcrop(xr1_0_m,[k kk windowz-1 windowz-1]);
                xr2_m=imcrop(xr2_0_m,[k kk windowz-1 windowz-1]);
                xr1=imcrop(xr1_0,[k kk windowz-1 windowz-1]);
                xr2=imcrop(xr2_0,[k kk windowz-1 windowz-1]);
                
                if saturation==1
                    xr1_m(xr1_m>si/100*255)=255;
                    xr2_m(xr1_m>si/100*255)=255;
                    
                    if j==1
                        imwrite(xr1_m,fullfile(imdir,['im_sat',num2str(si),'.png']));

                    end
                    
                end
                
                if gausblur==1
                    xr1_m=imgaussfilt(xr1_m,2);
                    xr2_m=imgaussfilt(xr2_m,2);
                end

                
                %convert uint8 variables to double precision variables so that
                %they can be used for calculations
                
                [X,Y]=meshgrid(1:size(xr1_m,2),1:size(xr1_m,1));
                
                if mbf==1
                    xr1d_m=double(xr1_m)-mean2(double(xr1_m));
                    xr2d_m=double(xr2_m)-mean2(double(xr2_m));
                    
                else
                    xr1d_m=double(xr1_m);
                    xr2d_m=double(xr2_m);
                end
                
                
                if dofilter==1
                    
                    %
                    %                     filter=exp(-((X-size(xr1,1)/2-1).^2./(2*(size(xr1,1)./6).^2)+...
                    %                         (Y-size(xr1,2)/2-1).^2./(2*(size(xr1,2)./6).^2)));
                    %                     x1f=filter.*xr1d;
                    %                     x2f=filter.*xr2d;
                    % %
                    filtersize=round(min(30,(windowz-2)/2));
                    PSF = fspecial('gaussian',30,5);
                    x1f = edgetaper(xr1d_m,PSF);
                    x2f = edgetaper(xr2d_m,PSF);
                else
                    x1f=xr1d_m;
                    x2f=xr2d_m;
                end
                
                if maf==1
                    x1f=x1f-mean2(x1f);
                    x2f=x2f-mean2(x2f);
                end
                
                
                F_CC=zeros(size(xr1_m));
                G1=fft2(x1f);
                G2=fft2(x2f);
                
                F_CC=G1.*conj(G2);
%                 F_AC=sqrt(G1.*conj(G1).*G2.*conj(G2));
                %%G1.*conj(G1);
                F_AC=0.5*(G1.*conj(G1)+G2.*conj(G2));%
                F_GCC=F_CC./abs(F_CC);
                F_GCC_adib=F_CC./F_AC;
                
                
                if ensemf==1
                    SCC_ensembleFrr=SCC_ensembleFrr+(F_CC);
                    AC_ensembleFrr=AC_ensembleFrr+(F_AC);
                    l=l+1;
                    
                    if abs(F_GCC)>0
                        PDF_GCCrr=PDF_GCCrr+F_GCC;
                        gcc_counter=gcc_counter+1;
                    end
                    if abs(F_GCC_adib)>0
                        PDF_Adibrr=PDF_Adibrr+F_GCC_adib;
                        adib_counter=adib_counter+1;
                    end
                    
                else
                    CC=fftshift(ifft2(F_CC));
                    AC=fftshift(ifft2(F_AC));
                    GCC=fftshift(ifft2(F_GCC));
                    GCC_adib=fftshift(ifft2(F_GCC_adib));
                    l=l+1;
                    SCC_ensembler=SCC_ensembler+CC;
                    AC_ensembler=AC_ensembler+AC;
                    if abs(GCC)>0
                        PDF_GCCr=PDF_GCCr+abs(GCC);
                        gcc_counter=gcc_counter+1;
                    end
                    if abs(GCC_adib)>0
                        PDF_Adibr=PDF_Adibr+abs(GCC_adib);
                        adib_counter=adib_counter+1;
                    end
                end
                
                if ensemf==1
                    ensgrids=ensgrids+1;
                    GCC_ensemble=(fftshift(ifft2(SCC_ensembleFrr./abs(SCC_ensembleFrr))));
                    SCC_ensembler=fftshift(ifft2(SCC_ensembleFrr));
                    AC_ensembler=fftshift(ifft2(AC_ensembleFrr));
                    PDF_GCCr=fftshift(ifft2(PDF_GCCrr));
                    PDF_Adibr=fftshift(ifft2(PDF_Adibrr));
                    pdfofensemble_Adib=(fftshift(ifft2(fft2(SCC_ensembler)./fft2(AC_ensembler))));
                    meannoiselevel(CorCounter)=mean(abs(pdfofensemble_Adib(:)));
                    integrallevel(CorCounter)=trapz(pdfofensemble_Adib(:));
                    L_infinity(ensgrids)=norm(pdfofensemble_Adib, inf);
                    if L_infinity(ensgrids)==inf || isnan(L_infinity(ensgrids)) 
                        try
                        L_infinity(ensgrids)=L_infinity(ensgrids-1);
                        catch
                        L_infinity(ensgrids)=1;
                        end
                    end
                    
%                     if convergencecriteria==1 && j>50+imstart && converged==0 %|| j==imstart+N-i-1
                        if convergencecriteria==1 && ensgrids>50 && converged==0 
%                         L_infinity(ensgrids);
                        convx=(1:ensgrids); % Explanatory variable
                        f = @(F,x) F(1)./sqrt(x)+F(2);
                        lb = [0,0];
                        ub = [inf,inf];
                        beta0 = [10,10];
                        opts = optimset('Display','off');
                        [FFF,resnorm] =lsqcurvefit(f,beta0,convx,L_infinity,lb,ub,opts);
                        eeerrr(ensgrids)=abs(L_infinity(ensgrids)-FFF(2))/FFF(2)*100;
                        FFFF(ensgrids)=FFF(2);
%                         FFF(2)

                        if  mean(eeerrr(ensgrids-49:ensgrids))<10 && std(FFFF(ensgrids-49:ensgrids))<0.1*mean(FFFF(ensgrids-49:ensgrids)) 
                            converged=CorCounter
                            figure
                            y_predicted=f(FFF,convx);
                            plot(1:ensgrids,L_infinity(1:ensgrids),'*');
                            hold on
                            plot(convx,y_predicted,'g','LineWidth',3);
                            plot([ensgrids ensgrids],[0 max(max(L_infinity),max(y_predicted))],'k','LineWidth',3)
                            axis tight
                            legend('data','fit','convergence')
%                             set(gca,'XScale','log')
                            pltname=strcat(resave,'/timelag',num2str(i),'convergence',num2str(CorCounter),'.fig')
                            saveas(gcf,pltname)
                            break
                        end
                        if  converged
                            breakflag=1
                            break
                        end
                    end
                    
                    
                else
                    GCC_ensemble=(ifft2(fft2(SCC_ensembler)./abs(fft2(SCC_ensembler))));
                    pdfofensemble_Adib=(fftshift(ifft2(fft2(SCC_ensembler)./fft2(AC_ensembler))));   %%%PDF using iPED
                    meannoiselevel(CorCounter)=mean(abs(pdfofensemble_Adib(:)));
                    integrallevel(CorCounter)=trapz(pdfofensemble_Adib(:));
                    L_infinity(CorCounter)=norm(pdfofensemble_Adib, inf);
                end
                
                xr1d=double(xr1)-mean2(double(xr1));
                xr2d=double(xr2)-mean2(double(xr2));
                F_CC_kahler=zeros(size(xr1_m));   %no filter for kahler+ ensemble in real
                G1=fft2(xr1d);
                G2=fft2(xr2d);
                F_CC_kahler=G1.*conj(G2);
                F_AC_kahler=(G1.*conj(G1));%0.5*(G1.*conj(G1)+G2.*conj(G2));%G1.*conj(G1);
                CC_kahler=fftshift(ifft2(F_CC_kahler));
                AC_kahler=fftshift(ifft2(F_AC_kahler));
                SCC_ensemble_kahler=SCC_ensemble_kahler+CC_kahler;
                AC_ensemble_kahler=AC_ensemble_kahler+AC_kahler;
            end
            if  converged
                breakflag=1
                break
            end
        end

        
        
        if Montecarlo==1
            SCC_ensemble = SCC_ensembler;%./(j);
            AC_ensemble = AC_ensembler;
            PDF_Adib=(PDF_Adibr)./(adib_counter);
            pdfofensemble_Adib=(fftshift(ifft2(fft2(SCC_ensemble)./fft2(AC_ensemble))));
            
            try
                [Diffusion_pdfofensemble_Adib(CorCounter)]= difes(pdfofensemble_Adib);
            catch
                [Diffusion_pdfofensemble_Adib(CorCounter)]=0;
            end
            try
                [Diffusion_PPDF_GCC(CorCounter)]= difes(PDF_GCCr);
            catch
                [Diffusion_PPDF_GCC(CorCounter)]= 0;
            end
            %                     try
            %                         [Diffusion_PDF_corrected(j)]= difes(PDF_corrected);
            %                     catch
            %                         [Diffusion_PDF_corrected(j)]= 0;
            %                     end
            %
            try
                [Diffusion_PDF_Adib(CorCounter)]= difes(PDF_Adibr);
            catch
                [Diffusion_PDF_Adib(CorCounter)]= 0;
            end
            [Dcc]= difes( SCC_ensemble_kahler);
            [Dac]= difes( AC_ensemble_kahler);
            Diffusion_kahler(CorCounter)=Dcc-Dac;
            if rem(CorCounter,100)==0
                savecor  = strcat(resave,'/convergecheck',num2str(N),'_',num2str(i),'_',num2str(CorCounter),'window',num2str(windowz),'_sat',num2str((si)));
                save(savecor,'Diffusion_pdfofensemble_Adib','Diffusion_kahler','Diffusion_PPDF_GCC','Diffusion_PDF_Adib')
            end
%             if getappdata(h,'canceling')
%                 breakflag=1;
%                 break
%                 
%             end
        end
        if  converged
            breakflag=1
            break
        end
%         waitbar(j/(N),h,strcat('Image ', {' '}, num2str(j),{' '}, 'out of ',{' '},num2str(N-i),{' '},' .',{' '},...
%             'Time lag', {' '},num2str(i),' out of',{' '},num2str(tmax), ' is being done. '))%{sprintf('\n')}
    end
    
    %     savecor  = strcat(resave,'\convergecheck',num2str(N),'_',num2str(i),'_',num2str(j),'_sat',num2str((si)));
    %     save(savecor,'pid','kahler')
    savemeannoise = strcat(resave,'/meannoise',num2str(CorCounter),'_sat',num2str((si)),'window',num2str(windowz),'_',num2str(i));
    save(savemeannoise,'meannoiselevel','integrallevel','L_infinity')
    
    
    if ensemf==1
        GCC_ensemble=abs(fftshift(ifft2(SCC_ensembleFrr./abs(SCC_ensembleFrr))));
        SCC_ensembler=fftshift(ifft2(SCC_ensembleFrr));
        AC_ensembler=fftshift(ifft2(AC_ensembleFrr));
        PDF_GCCr=fftshift(ifft2(PDF_GCCrr));
        PDF_Adibr=fftshift(ifft2(PDF_Adibrr));
        
    else
        GCC_ensemble=abs(ifft2(fft2(SCC_ensembler)./abs(fft2(SCC_ensembler))));
    end
    
    SCC_ensemble = SCC_ensembler;%;./l;%/(N-i);
    AC_ensemble = AC_ensembler;%./l;%/(N-i);
    PDF_GCC = PDF_GCCr./(gcc_counter);
    PDF_corrected=PDF_GCCr./(filter.^2);
    PDF_Adib=(PDF_Adibr./(adib_counter));
    pdf_iPED=(fftshift(ifft2(fft2(SCC_ensemble)./fft2(AC_ensemble))));
    

    
    savecor  = strcat(resave,'/pdf_iPED_',num2str(CorCounter),'window',num2str(windowz),'_','timelag',num2str(i));
    save(savecor,'pdf_iPED')

    
    
    if fittt==1 && Montecarlo==0
        
        try
            [Diffusion_pdfofensemble_Adib(i)]= difes(pdfofensemble_Adib);
        catch
            [Diffusion_pdfofensemble_Adib(i)]=0;
        end
        try
            [Diffusion_PPDF_GCC(i)]= difes(PDF_GCCr);
        catch
            [Diffusion_PPDF_GCC(i)]= 0;
        end
        try
            [Diffusion_PDF_corrected(i)]= difes(PDF_corrected);
        catch
            [Diffusion_PDF_corrected(i)]= 0;
        end
        
        try
            [Diffusion_PDF_Adib(i)]= difes(PDF_Adibr);
        catch
            [Diffusion_PDF_Adib(i)]= 0;
        end
        [Dcc,dd1,dd2,xcenter,ycenter,amp]= difes(SCC_ensemble_kahler);
        [Dac,Dacx,Dacy]= difes(AC_ensemble_kahler);
        Diffusion_kahler(i)=Dcc-Dac;
    end
    
    if Montecarlo==1
        savecor  = strcat(resave,'/convergence');
        save(savecor)
    end
end


if fittt==0
    diffusionestimate.pdfofensemble_Adib=0;
    delete(h)
elseif Montecarlo==0
    diffusionestimate.kahler=Diffusion_kahler;
    diffusionestimate.kahler_amp=amp;
    diffusionestimate.kahler_x=xcenter;
    diffusionestimate.kahler_y=ycenter;
    diffusionestimate.pdfofensemble_Adib=Diffusion_pdfofensemble_Adib;
    diffusionestimate.PDF_Adib= Diffusion_PDF_Adib;
    diffusionestimate.GCC=Diffusion_PPDF_GCC;
    diffusionestimate.GCC_corrected=Diffusion_PDF_corrected;
    diffusionestimate.particleACsize=4*sqrt(Dac);
    diffusionestimate.particleACdifX=(Dacx);
    diffusionestimate.particleACdifY=(Dacy);
    
    
    savecor  = strcat(resave,'/Diffestimate',num2str(CorCounter),'_',num2str(i),'window',num2str(windowz));
    save(savecor,'diffusionestimate')
    
    
%     delete(h)
end

toc
end





function [Diffusion_PDF,Diffusion_PDFx,Diffusion_PDFy,xcenter,ycenter,amp]= difes(Z)
% Z=abs(Z);
% Z=Z-min(Z(:));
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
%     F = @(xx,yy) x(1)*exp(-((xx-x(2)).^2/(2*x(3)^2) + (yy-x(4)).^2/(2*x(5)^2)));
ycenter=x(4);
xcenter=x(2);
pdfw=1/2*(x(3)+x(5));
Diffusion_PDF=pdfw^2/2;
Diffusion_PDFx=x(3)^2/2;
Diffusion_PDFy=x(5)^2/2;
amp=x(1);
end



