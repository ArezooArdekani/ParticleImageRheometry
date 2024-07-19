clc
% clear
close all


[Dpdf1,dpdf1x,dpdf1y,xc,yc,amp,F,x,Dpdf2]= ...
    difes(AC_ensemble);

particlesize=Dpdf1/sqrt(2)
% figure;imagesc(pdfofensemble_Adib);
% Diffusion=Dpdf
% [Dcc,dd1,dd2,xcenter,ycenter,amp]= difes(SCC_ensemble)
% [Dac,Dacx,Dacy,xcenter,ycenter,amp]= difes(AC_ensemble)
% Diffusion_kahler=Dcc-Dac



function [Diffusion_PDF,Diffusion_PDFx,Diffusion_PDFy,xcenter,...
    ycenter,amp,F,x,Diffusion_PDF2]= difes(Z)
Z=abs(Z);
Z=Z-min(Z(:));
figure(10)
plot(Z)

[X,Y]=meshgrid(1:size(Z,2),1:size(Z,1));
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
MdataSize = min(size(X,1),size(Y,2));
lb = [0,size(Z,2)/2+1,0,size(Z,2)/2+1,0,0,0,0,0,0,...
      0,size(Z,2)/2+1,0,size(Z,2)/2+1,0,0,0,0,0,0,...
      0,size(Z,2)/2+1,0,size(Z,2)/2+1,0,0,0,0,0,0,...
      0,size(Z,2)/2+1,0,size(Z,2)/2+1,0,0,0,0,0,0];
ub = [realmax('double'),size(Z,2)/2+1,(MdataSize)^2,MdataSize,(MdataSize)^2,realmax('double'),MdataSize,(MdataSize)^2,MdataSize,(MdataSize)^2, ...
    realmax('double'),size(Z,2)/2+1,(MdataSize)^2,MdataSize,(MdataSize)^2,realmax('double'),MdataSize,(MdataSize)^2,MdataSize,(MdataSize)^2, ...
    realmax('double'),size(Z,2)/2+1,(MdataSize)^2,MdataSize,(MdataSize)^2,realmax('double'),MdataSize,(MdataSize)^2,MdataSize,(MdataSize)^2, ...
    realmax('double'),size(Z,2)/2+1,(MdataSize)^2,MdataSize,(MdataSize)^2,realmax('double'),MdataSize,(MdataSize)^2,MdataSize,(MdataSize)^2];
x0 = [max(Z(:)),size(Z,2)/2+1,5,size(Z,1)/2+1,5,max(Z(:)),size(Z,2)/2+1,1,size(Z,1)/2+1,1,...
    max(Z(:)),size(Z,2)/2,5,size(Z,1)/2,5,max(Z(:)),size(Z,2)/2,5,size(Z,1)/2,5,...
    max(Z(:)),size(Z,2)/2,5,size(Z,1)/2,5,max(Z(:)),size(Z,2)/2,5,size(Z,1)/2,5,...
    max(Z(:)),size(Z,2)/2,55,size(Z,1)/2,50,max(Z(:)),size(Z,2)/2,55,size(Z,1)/2,50];
% opts = optimset('TolX',1e-30);
opts = optimoptions('lsqcurvefit','FunctionTolerance',1.0000e-12);%,'Algorithm','levenberg-marquardt');
[x,resnorm,residual,exitflag] = lsqcurvefit(@D5GaussFunction,x0,xdata,Z,lb,ub,opts);
F = @(xx,yy) x(1)*exp(   -((xx-x(2)).^2/(2*x(3)^2) + (yy-x(4)).^2/(2*x(5)^2) )    )+...
    x(6)*exp(   -((xx-x(7)).^2/(2*x(8)^2) + (yy-x(9)).^2/(2*x(10)^2) )    )+...
    x(11)*exp(   -((xx-x(12)).^2/(2*x(13)^2) + (yy-x(14)).^2/(2*x(15)^2) )    )+...     
    x(16)*exp(   -((xx-x(17)).^2/(2*x(18)^2) + (yy-x(19)).^2/(2*x(20)^2) )    );%+
ycenter=x(4)
xcenter=x(2)
pdfw=1/2*(x(3)+x(5));
Diffusion_PDF=pdfw^2/2;
Diffusion_PDFx=x(3)^2/2;
Diffusion_PDFy=x(5)^2/2;
size1=4*pdfw/sqrt(2)

amp=x(1)
amp2=x(6)
ycenter2=x(7)
xcenter2=x(9)
pdfw2=1/2*(x(8)+x(10));
Diffusion_PDF2=pdfw2^2/2;
size2=4*pdfw2/sqrt(2)
amp3=x(11)
ycenter3=x(12)
xcenter3=x(14)
pdfw3=1/2*(x(13)+x(15));
Diffusion_PDF3=pdfw3^2/2;
size3=4*pdfw3/sqrt(2)

% amp4=x(16)
% ycenter4=x(17)
% xcenter4=x(19)
pdfw4=1/2*(x(18)+x(20));
% Diffusion_PDF4=pdfw4^2/2
size4=4*pdfw4/sqrt(2)

aaaa=F(X,Y);
% figure;imagesc(aaaa);


h = figure();
subplot(4,4, [4])
imagesc(aaaa)
hold on
axis tight
subplot(4,4, [5,6,7,9,10,11,13,14,15])
imagesc(X(1,:),Y(:,1)',Z)
% axis equal
% xlim([size(Z,2)/2-30 size(Z,2)/2+30])
% ylim([size(Z,2)/2-30 size(Z,2)/2+30])

set(gca,'YDir','reverse')
colormap('jet')
hold on
[row, col] = find(ismember(Z, max(Z(:))));
m =0;% Point slope formula
b = (-m*x(2) + row);
xvh = 1:size(Z,2);
yvh = xvh*m + b;
hPoints = Z(col,:);%size(Z,1)/2+1,:);
% generate points along vertical axis
mrot = -m;
brot = (mrot*x(4) - col);
yvv = 1:size(Z,1);
xvv = yvv*mrot - brot;
vPoints = Z(:,row);%size(Z,2)/2+1);%interp2(X,Y,Z,xvv,yvv);
% plot lins
plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r')
plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'g')
xdatafit = linspace(1 ,MdataSize,10000);
hdatafit = F(xdatafit,col);%x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
vdatafit = F(row,xdatafit);% x(1)*exp(-(xdatafit-x(4)).^2/(2*x(5)^2));
subplot(4,4, [1:3])
xposh = (xvh-x(2))/1+x(2);% correct for the longer diagonal if fi~=0
% plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black')
% hold on
plot(xposh,Z')
hold on
plot(xdatafit,hdatafit,'black','LineWidth',3)
set(gca,'YScale','log')

axis tight
% xlim([size(Z,2)/2-30 size(Z,2)/2+30])

xposv = (yvv-x(4))/cos(0)+x(4);
subplot(4,4,[8,12,16])
plot(vPoints,xposv,'g.',vdatafit,xdatafit,'black')
set(gca,'YScale','log')

axis tight
% ylim([size(Z,2)/2-30 size(Z,2)/2+30])
end
