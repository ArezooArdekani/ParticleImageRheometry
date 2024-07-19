function [D3p]=fit2d_3p(CCm,tr,dimen,Diff)

% [row, col] = find(ismember(CCm, max(CCm(:))))
% row=129;
% col=129;
[row, col] = find(ismember(CCm, max(max(CCm(120:130,120:130)))))
% % CCm=CCm(:,col)-mean2(CCm(:,1));
% xx=1:size(CCm,1);
% cross=CCm(:,col);%-mean(CCm(1:100,col));
% figure(1)
% plot(cross)
% hold on
% 
% % cross_trunc=[cross(1:row-3);cross(row+3:end)];
% % xcross_trun=[x(1:row-3) x(row+3:end)];
% 
% x=1:size(CCm,1);
% cross_trunc=CCm(row-tr:row+tr,col)-min(CCm(row-tr:row+tr,col));
% xcross_trun=x(row-tr:row+tr);
% 
% figure(2)
% plot(xcross_trun,cross_trunc,'bs')
% hold on
% % plot(x,cross,'r.')
% 
% % f= fit(xx',xcross,'gauss1')
% f= fit(xcross_trun',cross_trunc,'gauss2');
% Ax=coeffvalues(f);
% cx(1)=Ax(3)/sqrt(2);
% Amp(1)=Ax(1);
% try
%     cx(2)=Ax(6)/sqrt(2);
%     Amp(2)=Ax(4);
%     cx(3)=Ax(9)/sqrt(2);
%     Amp(3)=Ax(7)
%     cx(4)=Ax(12)/sqrt(2);
%     Amp(4)=Ax(10)
% end
% [M,I] = max(Amp);
% cx1=cx(I);
% CCm=CCm-min(CCm(:));
% CCm=CCm./(CCm(row,col));

A=[1 0 0 0 0;
   1 -1 0 1 0;
   1 1 0 1 0;
   1 0 -1 0 1;
   1 0 1 0 1];

Z=[ log(CCm(row,col));
    log(CCm(row,col-1));
    log(CCm(row,col+1));
    log(CCm(row-1,col));
    log(CCm(row+1,col))];

c=A^-1*Z;
% 
 syms f(x,y)
 fz =@(x,y) exp(c(1)+c(2).*y+c(3).*x+c(4).*y.^2+c(5).*x.^2);
for x=1:256
    for y=1:256
        z(x,y)=fz(x-col,y-row);
    end
end
c
sx=sqrt(1/abs(2*c(4)));
sy=sqrt(1/abs(2*c(5)));
D3p=(0.5*(sy+sx))^2/2;


% x=1:size(CCm,1);
% xcross=z(x,row);
% PDF_truncated =[x' xcross];
% f= fit(x',xcross,'gauss1');
% Ax=coeffvalues(f);
% cx1=Ax(3)/sqrt(2);
% try
%     cx2=Ax(6)/sqrt(2);
%     cx3=Ax(9)/sqrt(2);
%     cx4=Ax(12)/sqrt(2);
% end
figure(1)
xx=1:size(CCm,1);
cross=CCm(:,col);
plot(xx,z(xx,col),xx,cross,'*')

hold on
% plot(xx,z(xx,col))
% D=cx1^2/2




% figure(5)
% plot(f,x,cross)
% hold on
% figure
% plot(f,x,CCm(x,row))
% hold on
% figure(5)
% plot(f,xcross_trun,cross_trunc)
% hold on


if dimen==1


y=1:size(CCm,1);
ycross=CCm(row,:)-mean(CCm(row,1:100));
% figure(3)
% plot(ycross,'bs')
ycross_trunc=[ycross(1:col-3) ycross(col+3:end)];
ycross_trun=[y(1:col-3) y(col+3:end)];
figure(3)
plot(ycross_trun',ycross_trunc','o');
hold on
plot(y,ycross,'r.')

% f= fit(xx',xcross,'gauss1')
f= fit(ycross_trun',ycross_trunc','gauss1');
Ay=coeffvalues(f);
cy(1)=Ay(3)/sqrt(2);
Amp(1)=Ay(1);
try
    cy(2)=Ay(6)/sqrt(2);
    Amp(2)=Ay(4);
    cy(3)=Ay(9)/sqrt(2);
    Amp(3)=Ay(7)
    cy(4)=Ay(12)/sqrt(2);
    Amp(4)=Ay(10)
end
[M,I] = max(Amp);
cy1=cy(I);

% figure(5)
% plot(f,x,cross)
% hold on
% figure
% plot(f,x,CCm(x,row))
% hold on
figure(6)
plot(f,ycross_trun,cross_trunc)
hold on


else 
    cy1=0;

end









% A=[1 0 0 0 0;
%    1 -1 0 1 0;
%    1 1 0 1 0;
%    1 0 -1 0 1;
%    1 0 1 0 1];
% CCm=CCm-mean2(CCm(1:10,col));
% Z=[ log(CCm(row,col));
%     log(CCm(row,col-1));
%     log(CCm(row,col+1));
%     log(CCm(row-1,col));
%     log(CCm(row+1,col))];
% 
% c=A^-1*Z;
% 
% % syms f(x,y)
% fz =@(x,y) exp(c(1)+c(2).*y+c(3).*x+c(4).*y.^2+c(5).*x.^2);
% 
% for x=1:256
%     for y=1:256
%         z(x,y)=fz(x-col,y-row);
%     end
% end
% 
% 
% x=1:size(CCm,1);
% xcross=z(x,row);
% PDF_truncated =[x' xcross];
% f= fit(x',xcross,'gauss1');
% Ax=coeffvalues(f);
% cx13p=Ax(3)/sqrt(2);
% try
%     cx23p=Ax(6)/sqrt(2);
%     cx33p=Ax(9)/sqrt(2);
%     cx43p=Ax(12)/sqrt(2);
% end
% % plot(f,'g',x,xcross,'gs')
% D=cx1^2/2;

