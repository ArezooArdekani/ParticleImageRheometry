function [cx1,cy1]=fit2d_filtered(CCm,tr)

% [row, col] = find(ismember(CCm, max(CCm(:))))
row=129
col=129
% CCm=CCm(:,col)-mean2(CCm(:,1));
x=1:size(CCm,1);
xcross=CCm(row-tr:row+tr,col)-min(CCm(row-tr:row+tr,col));
xx=x(row-tr:row+tr);

% plot(xx,xcross,'bs')
% hold on
% plot(xx,xcross,'r.')

% f= fit(xx',xcross,'gauss1')
f= fit(xx',xcross,'gauss1');
Ax=coeffvalues(f);
cx(1)=Ax(3)/sqrt(2);
Amp(1)=Ax(1);
try
    cx(2)=Ax(6)/sqrt(2);
    Amp(2)=Ax(4);
    cx(3)=Ax(9)/sqrt(2);
    Amp(3)=Ax(7);
    cx(4)=Ax(12)/sqrt(2);
    Amp(4)=Ax(10);
end
cx
[M,I] = max(Amp);
cx1=cx(I);

% figure(5)
% plot(f,x,cross)
% hold on
% figure
% plot(f,x,CCm(x,row))
% hold on
figure(5)
plot(f,xx,xcross)
hold on





y=1:size(CCm,1);
ycross=CCm(row,row-tr:row+tr)-min(CCm(row,row-tr:row+tr));
% figure(3)
% plot(ycross,'bs')
yy=y(row-tr:row+tr);

% figure(3)
% plot(yy',ycross','o');
% hold on
% plot(yy,ycross,'r.')

% f= fit(xx',xcross,'gauss1')
f= fit(yy',ycross','gauss1');
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
cy
[M,I] = max(Amp);
cy1=cy(I);

% figure(5)
% plot(f,x,cross)
% hold on
% figure
% plot(f,x,CCm(x,row))
% hold on
figure(6)
plot(f,yy,ycross)
hold on


figure(10)
mesh(CCm(row-tr:row+tr,col-tr:col+tr))
















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

