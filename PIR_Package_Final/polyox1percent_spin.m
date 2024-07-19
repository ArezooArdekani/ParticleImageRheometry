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
for tracking=0:1:1
pir=[];
pird=[];
pirmsd=[];
comp_V=[];
 for casesi=1:size(caseha,2)
     cases=caseha(casesi);
 if tracking==1   
load(strcat('P:\Projects\uPIV\polyox4m_1percent_spinningdisc_20x_500nm\msd_of_',num2str(cases),'microrheology_1P.mat'));
ms2=msd';
clear msd
scale=6.5/20*10^-6;
 else
resave=strcat('P:\Projects\uPIV\RESULTS_polyox4m_1percent_spinningdisc_20x_500nm\r',num2str(cases),'\results_1999');
%W:\Projects\Microscale_measurments\micro-rheology\confocal\polyox\low_conc\r1\results_full';
% resave='W:\Projects\Microscale_measurments\micro-rheology\confocal\polyox\high_conc\r1\results_full';
%'E:\BME_Confocal\Confocal_data\100nm_3k\gauss1_backsub\results'
% filetype='.mat'
% folderPath1 = resave;
% txtpatternpdf = fullfile(folderPath1, filetype);
% dinfoAC = dir(txtpatternpdf)
% 
% [~, reindex] = sort( str2double( regexp( {dinfoAC.name}, '\d+', 'match', 'once' )))
% dinfoAC = dinfoAC(reindex) ;
% 
% pdf1=load(strcat(dinfoAC(1).folder,'/',dinfoAC(1).name));
    
filename='pdfofensemble_Adib_999_sat100window64_'%'PDF_ensemble_5000_'%
% % % 
try
% resave=strcat('W:\Projects\Microscale_measurments\micro-rheology\pir_new\polyox1percent_200nm_23C\r',num2str(cases),'\results_999');
load(strcat(resave,'\msd.mat'))
catch
[ms2]=PDF2MSD(adiboutlier,matout,resave,filename,tmin,tmax,2000)

end
scale=3.5/20*10^-6;
 end

ms2=ms2*(scale)^2;%(0.325*10^-6)^2;%(0.433)^2*10^-12;%(0.16)^2*10^-12;%*(10/60)^2*10^-12;% *0.11^2*10^-12% %*(10/60)^2*10^-12%*0.11^2*10^-12%(10/60)^2*10^-12;  %[m^2]
dt=1/9.335;%7.603%10.3%9.66;%33.76%13.517%33.76%13.517%17.1%9.1%49.319%9.1;%3158%2.002%0.08049%2.002%0.08049%1/26.665%1./17.172;%0.029%%0.029%1/26.665%23.043%30%0.08049%2; 2%  %[S]

r=0.5*500*10^-9  %196*10^-12%0.5*500*10^-9  [m]
timelags=(tmin:tmax).*dt;
msd=ms2(tmin:tmax);%4*2.16184665048e-12*timelags%



% % handles.pltw=plot(timelags,msd,'r*');
% % grid on
% % xlabel('time lag');
% % ylabel('MSD');
% figure(2)
% plot(timelags,msd,'r*');
% grid on
% xlabel('time lag [s]');
% ylabel('MSD [m^2]');













Kb=1.38064852e-23;% m2 kg s-2 K-1
T=22.5 +273.15; %K
% mu= 10e-4; %pa.s
% a=100e-9; %m
% % D=Kb*T/(6*pi*mu*r); %[m^2./s] 
% % Dt=2;
% r=0.5*500*10^-9  
%Kb*T/(6*pi*mu*Dt); %[m^2./s]
% % Dt=Kb*T/(6*pi*mu*a)
% % % Dum2=D*10^12
% % % scale=(1/6)*10^-6; % m/px
% % % dt=10; %frme/s
% % % D_px_fr=(D/(scale^2))*1./dt

% % % dmea=0.52*(scale^2)*dt
% % % plot(,)
% timelags=1:1000%q(:,1)'
% msd=6*Dt*timelags.^1% q(:,2)';
% % % % % msd(101:999)=msd(101:999)+4*Dt*timelags(101:999).^0.5
%%
figure(10)
subplot(1,2,1)
plot(timelags,msd,'b>','LineWidth',3)
xlabel('\tau [s]')
ylabel('MSD')
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
timelagsr=timelags;
msdr=msd;
% [timelagsr,msdr]=rbfmsd(timelags,msd)
% timelagsr = linspace(min(timelags),max(timelags),50);
% msdr = spline(timelags,msd,timelagsr);


% l_msd=log(msdr);
% l_t=log(timelagsr);
% d_msd_dt=gradient(l_msd)./gradient(l_t);
% 
% 
% xd = diff([l_t(3),l_t,l_t(end-2)]);
% ud = diff([l_msd(3),l_msd,l_msd(end-2)]);
% dudx = (ud(1:end-1)./xd(1:end-1).*xd(2:end) ...
%           + ud(2:end)./xd(2:end).*xd(1:end-1)) ...
%           ./ (xd(2:end)+xd(1:end-1));
% d_msd_dt=dudx;

% w=2*pi./timelagsr;
% omega=w;
% w=[omega 1];
% for i=1:size(w,2)-1
%     i
% G_w(i)=Kb*T/(pi*r*msdr(i)*gamma(1+d_msd_dt(i)));
% G_p(i)=abs(G_w(i))*cos(pi*d_msd_dt(i)/2);
% G_dp(i)=abs(G_w(i))*sin(pi*d_msd_dt(i)/2);
% n_c(i)=sqrt((abs(G_p(i)./(w(i)))^2+abs(G_dp(i)./(w(i)))^2));%./(w(i)^2));
% end
if tracking ==1
[msdtau] = making_logarithmically_spaced_msd_vs_tau(msdr*10^12,timelags,max(timelags))
[ omega,Gs,Gp,Gpp, dd, dda ] = calc_G(msdtau(:,1),msdtau(:,2),r*10^6,2,T,0.03,2)%0.3,1.5)
else
[msdtau] = making_logarithmically_spaced_msd_vs_tau(msdr*10^12,timelags,max(timelags))
[ omega,Gs,Gp,Gpp, dd, dda ] = calc_G(msdtau(:,1),msdtau(:,2),r*10^6,2,T,0.03,2)%0.3,1.5)   
end
G_w=Gs;%Kb*T/(pi*r*msdr(i)*gamma(1+d_msd_dt(i)));
G_p=Gp;%abs(G_w(i))*cos(pi*d_msd_dt(i)/2);
G_dp=Gpp;%abs(G_w(i))*sin(pi*d_msd_dt(i)/2);
n_c=Gs;%sqrt((abs(G_p(i)./(w(i)))^2+abs(G_dp(i)./(w(i)))^2));%./(w(i)^2));
w=[omega 1];

% close all
%%
figure(1)
subplot(1,2,1)
plot(timelags,msd,'b>','LineWidth',3)
hold on;
% plot(timelagsr,msdr,'k','LineWidth',1)
xlabel('\tau [s]')
ylabel('MSD [m^2]')
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
legend({'MSD','RBF representation'})


figure(1)
subplot(1,2,2)
% plot(w(1:end-1),G_w,'s','LineWidth',3)
hold on
% yyaxis left
% plot(w(1:end-1),n_c,'ko-','LineWidth',3)
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
ylabel('Complex viscosity [Pa.s]')
% % ylim([10^-2 30])
% yyaxis right
if tracking ==1
plot(w(1:end-1),G_p,'ms','LineWidth',2,'MarkerSize',10)
hold on
plot(w(1:end-1),G_dp,'go','LineWidth',2,'MarkerSize',10)
legend('Gp tracking','G" tracking')
xlabel('\omega [rad/s]')
else
plot(w(1:end-1),G_p,'ms','LineWidth',2,'MarkerSize',10)
hold on
plot(w(1:end-1),G_dp,'go','LineWidth',2,'MarkerSize',10)
legend('Gp PIR','G" PIR')
xlabel('\omega [1/s]')
end
% ww=flip(w(1:end-1));
% n_c=sqrt((abs(G_p./(ww)).^2+(abs(G_dp./(ww))).^2));
% plot(w(1:end-1),n_c,'r>','LineWidth',2)
% ylim([10^-2 30])
% % legend('G*','G prime','G doubleprime','complex viscosity')
% legend('\eta','G''','G"')

% ylabel('G'' and G" [Pa]')
% xlim([10^-2 20])
box on
% ylim([10^-4 10^-1])
% xlim([10^-1 10])
% xlim([1 30])
% xlim([0.9 10])

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
% set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
% legend({'MSD Tracking'})
else
errorbar(timelags,mean(pirmsd),std(pirmsd),'r<','LineWidth',3,'DisplayName','MSD PIR')
hold on;
xlabel('\tau [s]')
ylabel('MSD [m^2]')
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
% ({'MSD PIR'})
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
% legend('','G" Tracking')

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
% errorbar(w(1:end-1),mean(CV),std(CV),'k.','LineWidth',2,'MarkerSize',10,'DisplayName','Complex G PIR')
% legend('Gp PIR','G" PIR')
end

% ww=flip(w(1:end-1));
% n_c=sqrt((abs(G_p./(ww)).^2+(abs(G_dp./(ww))).^2));
% plot(w(1:end-1),n_c,'r>','LineWidth',2)
% ylim([10^-2 30])
% % legend('G*','G prime','G doubleprime','complex viscosity')
% legend('\eta','G''','G"')

xlabel('\omega [1/s]')
% ylabel('G'' and G" [Pa]')
% xlim([10^-2 20])
box on
% % ww=flip(w(1:end-1));
% % CV_PR=sqrt((abs(pir./(ww)).^2+(abs(pird./(ww))).^2));
% % errorbar(w(1:end-1),mean(CV_PR),std(CV_PR),'k>','LineWidth',2)
% legend('Gp','Gdp','complexviscosity')
% xlabel('\omega [1/s]')
% box on
set(gca,'FontSize',20,'LineWidth',3)
end

legend
% keyboard
%% TA data

%BME TA
freq=[0.1
    0.125892
    0.158489
    0.199525
    0.25119
    0.316225
    0.398109
    0.501187
    0.630956
    0.794337
    1
    1.25891
    1.58488
    1.99528
    2.51191
    3.16224
    3.98106
    5.01187
    6.30951
    7.94323
    10
    12.5893
    15.8489
    19.9528
    25.119
    31.6228
    39.811
    50.1187
    63.0956
    79.4333
    100
    ]

gprime=10^6*[[2.63988e-8
3.80610e-8
5.37179e-8
7.49195e-8
1.02574e-7
1.39616e-7
1.87728e-7
2.48578e-7
3.25322e-7
4.21428e-7
5.38384e-7
6.79133e-7
8.45508e-7
1.03732e-6
1.26185e-6
1.51648e-6
1.81249e-6
2.14246e-6
2.51866e-6
2.93499e-6
3.37027e-6
3.76292e-6
3.93624e-6
3.46963e-6
1.98213e-6
2.24889e-7
-8.89094e-7
4.05229e-6
2.99200e-5
8.83250e-5
3.52989e-4] [2.25716e-8
3.33939e-8
4.85898e-8
6.95714e-8
9.77535e-8
1.35247e-7
1.85106e-7
2.47928e-7
3.27454e-7
4.28042e-7
5.49423e-7
6.97460e-7
8.71614e-7
1.07435e-6
1.31067e-6
1.58175e-6
1.89100e-6
2.25050e-6
2.65780e-6
3.12945e-6
3.65270e-6
4.14864e-6
4.43988e-6
4.24673e-6
2.98637e-6
-1.41024e-7
1.77934e-6
8.68620e-6
4.09702e-5
1.18097e-4
2.73959e-4] [2.26892e-8
3.36348e-8
4.87611e-8
7.00045e-8
9.93918e-8
1.38248e-7
1.88743e-7
2.53980e-7
3.36038e-7
4.38747e-7
5.64121e-7
7.15088e-7
8.95496e-7
1.10370e-6
1.34882e-6
1.62944e-6
1.95729e-6
2.33021e-6
2.76221e-6
3.26588e-6
3.82201e-6
4.38358e-6
4.75507e-6
4.56138e-6
3.43251e-6
3.86313e-6
2.23250e-6
9.11905e-6
3.93207e-5
1.15945e-4
2.79258e-4]]

gdoubleprime=10^6*[[1.88994e-7
2.33204e-7
2.85995e-7
3.48853e-7
4.24340e-7
5.12132e-7
6.15074e-7
7.34370e-7
8.71019e-7
1.02740e-6
1.20470e-6
1.40548e-6
1.63028e-6
1.88094e-6
2.15789e-6
2.46335e-6
2.79802e-6
3.16801e-6
3.57238e-6
4.02327e-6
4.52389e-6
5.08467e-6
5.73610e-6
6.53755e-6
7.58357e-6
8.27362e-6
7.72063e-6
1.96547e-6
-1.49607e-5
-5.75790e-5
-2.71290e-4] [1.89038e-7
2.33974e-7
2.89048e-7
3.54456e-7
4.31799e-7
5.22826e-7
6.28972e-7
7.51833e-7
8.92780e-7
1.05259e-6
1.23539e-6
1.43993e-6
1.66862e-6
1.92473e-6
2.20506e-6
2.51545e-6
2.85427e-6
3.22684e-6
3.63754e-6
4.08375e-6
4.58900e-6
5.15469e-6
5.80638e-6
6.56717e-6
7.61712e-6
8.87778e-6
7.49575e-6
2.23251e-6
-1.63279e-5
-6.72174e-5
-4.05640e-4] [1.93074e-7
2.39340e-7
2.95282e-7
3.62236e-7
4.41477e-7
5.34520e-7
6.43141e-7
7.68007e-7
9.12045e-7
1.07517e-6
1.25889e-6
1.46722e-6
1.69933e-6
1.95792e-6
2.24152e-6
2.55339e-6
2.89556e-6
3.27166e-6
3.68178e-6
4.13358e-6
4.63221e-6
5.20736e-6
5.85312e-6
6.63056e-6
7.73994e-6
8.19424e-6
7.56469e-6
1.54741e-6
-1.57337e-5
-6.81698e-5
-3.74339e-4]]

figure(2)
subplot(1,3,2)
hold on
errorbar(freq,mean(gprime'),std(gprime'),'sr-','LineWidth',2,'DisplayName','Gp TA','MarkerFaceColor','r','MarkerSize',10)
xlabel('\omega [rad/s]')
ylabel('G [Pa]')
hold on
legend
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
ylim([10^-5 10^5])
subplot(1,3,3)
xlabel('\omega [rad/s]')
hold on
ylabel('G" [Pa]')
ylim([10^-5 10^5])
errorbar(freq,mean(gdoubleprime'),std(gdoubleprime'),'ob-','LineWidth',2,'DisplayName','G" TA','MarkerFaceColor','b','MarkerSize',10)
% errorbar(freq,mean(CV')',std(CV')','k>-','LineWidth',3,'DisplayName','Complex V_TA')
legend
set(gca,'FontSize',20,'LineWidth',3,'XScale','log','YScale','log')
% ylim([10^-5 10^5])