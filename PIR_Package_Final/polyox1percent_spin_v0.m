clc
clear
close all
adiboutlier=1;
matout=1;
tmin=1;
tmax=30%00;%50
caseha=[3 4 5]% 6 7 8]% 6]%6]%6 7 8]


 addpath 'W:\Projects\Microscale_measurments\micro-rheology\Parralel_RPM\microrheology-master\src\Moduli';
 addpath 'W:\Projects\Microscale_measurments\micro-rheology\Parralel_RPM\microrheology-master\src\MSD';

[ms2]=PDF2MSD(adiboutlier,matout,resave,filename,tmin,tmax,2000)

end
scale=3.5/20*10^-6;
 end

ms2=ms2*(scale)^2;%(0.325*10^-6)^2;%(0.433)^2*10^-12;%(0.16)^2*10^-12;%*(10/60)^2*10^-12;% *0.11^2*10^-12% %*(10/60)^2*10^-12%*0.11^2*10^-12%(10/60)^2*10^-12;  %[m^2]
dt=1/9.335;%7.603%10.3%9.66;%33.76%13.517%33.76%13.517%17.1%9.1%49.319%9.1;%3158%2.002%0.08049%2.002%0.08049%1/26.665%1./17.172;%0.029%%0.029%1/26.665%23.043%30%0.08049%2; 2%  %[S]

r=0.5*500*10^-9  %196*10^-12%0.5*500*10^-9  [m]
timelags=(tmin:tmax).*dt;
msd=ms2(tmin:tmax);%timelags%
Kb=1.38064852e-23;% m2 kg s-2 K-1
T=22.5 +273.15; %K
%%
figure(10)
subplot(1,2,1)
plot(timelags,msd,'b>','LineWidth',3)
xlabel('\tau [s]')
ylabel('MSD')
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
timelagsr=timelags;
msdr=msd;
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
pirmsd=[pirmsd;msd];
pir =[pir;G_p];
pird=[pird;G_dp];
comp_V=[comp_V;n_c];
end
%%
w=[omega 1];
figure(2)
subplot(1,3,1)
if tracking ==1
errorbar(timelags,mean(pirmsd),std(pirmsd),'b>','LineWidth',3,'DisplayName','MSD Tracking')
hold on;
xlabel('\tau [s]')
ylabel('MSD [m^2]')
else
errorbar(timelags,mean(pirmsd),std(pirmsd),'r<','LineWidth',3,'DisplayName','MSD PIR')
hold on;
xlabel('\tau [s]')
ylabel('MSD [m^2]')
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
end
legend

subplot(1,3,2)
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
ylabel('Complex viscosity [Pa.s]')
if tracking ==1
    loyolagreen = [0.8500, 0.3250, 0.0980];
        subplot(1,3,2)
    errorbar(w(1:end-1),mean(pir),std(pir),'s','Color', loyolagreen,'LineWidth',2,'MarkerSize',10,'DisplayName','Gp Tracking')
hold on
subplot(1,3,3)
hold on
errorbar(w(1:end-1),mean(pird),std(pird),'co','LineWidth',2,'MarkerSize',10,'DisplayName','G" Tracking')

else
subplot(1,3,2)
errorbar(w(1:end-1),mean(pir),std(pir),'ms','LineWidth',2,'MarkerSize',10,'DisplayName','Gp PIR')
hold on
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
subplot(1,3,3)
errorbar(w(1:end-1),mean(pird),std(pird),'go','LineWidth',2,'MarkerSize',10,'DisplayName','G" PIR')
hold on
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
CV=sqrt((abs(pir./(w(1:end-1))).^2+(abs(pird./(w(1:end-1)))).^2));
end

xlabel('\omega [1/s]')
box on
set(gca,'FontSize',20,'LineWidth',3)
end
