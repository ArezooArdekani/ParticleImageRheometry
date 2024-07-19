clc
% clear
close all
% 
% resave='E:\PID\phantom\194nm\dofilter1_gaussb0_mbf1maf0ensemf1';
% filename='pdfofensemble_Adib_200_';

for i=1:1
    i
    % filedir=strcat(resave,'\',filename,num2str(i),'_sat100.mat')
    % load(filedir);
    %     Z = AC_ensemble;%pdfofensemble_Adib;
    %     [row, col] = find(ismember(Z, max(Z(:))));
    %     x=[1:size(Z,1)/2 size(Z,1)/2+2:size(Z,1)];
    %     xcross=Z(x,row);
    %     f= fit(x',xcross,'gauss2');
    %     figure(100)
    %     plot(f,x,xcross)
    %     Ax=coeffvalues(f);
    %     cx1=Ax(3)/sqrt(2)
    %     cx2=Ax(5)/sqrt(2)
    %
    %     r1=f(size(Z,1)/2+1);
    %     xcross=Z(col,x);
    %     f= fit(x',xcross','gauss2');
    %     figure(2)
    %     plot(f,x,xcross)
    %     Ax=coeffvalues(f);
    %     cx1=Ax(3)/sqrt(2);
    %     r2=f(size(Z,1)/2+1);
    % %     pdfofensemble_Adib(size(Z,2)/2+1,size(Z,2)/2+1)=0.5*(r1+r2);
    %     AC_ensemble(size(Z,2)/2+1,size(Z,2)/2+1)=0.5*(r1+r2);
    
%     [Dpdf,dpdf1,dpdf2,xc,yc]= difes(pdfofensemble_Adib)
%     [Dpdf,dpdf1,dpdf2,xc,yc]= difes(PDF_GCC)
%     Diffusion=Dpdf;
%     [Dcc,dd1,dd2,xcenter,ycenter,amp]= difes(SCC_ensemble_kahler);
    [Dac,Dacx,Dacy,xcenter,ycenter,amp]= difes(AC_ensemble)
%     Diffusion_kahler=Dcc-Dac
%     
% eiped=100*abs(5-Dpdf)/Dpdf
% ekahler=100*abs(5-Diffusion_kahler)/Diffusion_kahler
% ipeds=Dpdf*0.11^2/0.065
% kahlers=Diffusion_kahler*0.11^2/0.065
end




function [Diffusion_PDF,Diffusion_PDFx,Diffusion_PDFy,xcenter,ycenter,amp]= difes(Z)
Z=abs(Z)-mean(Z(:));
% Z=Z;

% [row, col] = find(ismember(Z, max(Z(:))));
% x=[1:size(Z,1)/2 size(Z,1)/2+2:size(Z,1)];
% Z(size(Z,2)/2+1,size(Z,2)/2+1)=0.5*(Z(size(Z,2)/2-1,size(Z,2)/2-1)+Z(size(Z,2)/2,size(Z,2)/2));
% xcross=Z(x,row);
% f= fit(x',xcross,'gauss2');
% figure(100)
% plot(f,x,xcross)
% Ax=coeffvalues(f);
% cx1=Ax(3)/sqrt(2)
% cx2=Ax(5)/sqrt(2)
% 
% r1=f(size(Z,1)/2+1);
% xcross=Z(col,x);
% f= fit(x',xcross','gauss2');
% figure(2)
% plot(f,x,xcross)
% Ax=coeffvalues(f);
% cx1=Ax(3)/sqrt(2);
% r2=f(size(Z,1)/2+1);
% Z(size(Z,2)/2+1,size(Z,2)/2+1)=0.5*(r1+r2);


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
h = figure();
subplot(4,4, [5,6,7,9,10,11,13,14,15])
imagesc(X(1,:),Y(:,1)',Z)
% axis equal
xlim([size(Z,2)/2-30 size(Z,2)/2+30])
ylim([size(Z,2)/2-30 size(Z,2)/2+30])

set(gca,'YDir','reverse')
colormap('jet')
hold on
m =0;% Point slope formula
b = (-m*x(2) + x(4));
xvh = 1:size(Z,2);
yvh = xvh*m + b;
hPoints = Z(size(Z,1)/2+1,:);
% generate points along vertical axis
mrot = -m;
brot = (mrot*x(4) - x(2));
yvv = 1:size(Z,1);
xvv = yvv*mrot - brot;
vPoints = Z(:,size(Z,2)/2+1);%interp2(X,Y,Z,xvv,yvv);
% plot lins
plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r')
plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'g')
xdatafit = linspace(1 ,MdataSize,10000);
hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
vdatafit = x(1)*exp(-(xdatafit-x(4)).^2/(2*x(5)^2));
subplot(4,4, [1:3])
xposh = (xvh-x(2))/1+x(2);% correct for the longer diagonal if fi~=0
plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black')
xposv = (yvv-x(4))/cos(0)+x(4);
axis tight
xlim([size(Z,2)/2-30 size(Z,2)/2+30])

subplot(4,4,[8,12,16])
plot(vPoints,xposv,'g.',vdatafit,xdatafit,'black')
axis tight
ylim([size(Z,2)/2-30 size(Z,2)/2+30])
end
