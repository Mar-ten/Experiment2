function [ ] = contourCompare( data )
%contourCompare plots comparisons between theoretical and experimental
%displacement fields. To run, first execute Mode1Analysis.m then use
%contourCompare(data) where data can be data1-data8 
data.x=data.x*1E3;
data.y=data.y*1E3;
data.u=data.u*1E3;
data.v=data.v*1E3;
data.uT=data.uT*1E3;
data.vT=data.vT*1E3;


figure;
hold on
contour(data.y,data.x,data.u,'b','ShowText','on')
contour(data.y,data.x,data.uT,'.-k','ShowText','on')
plot(0,0,'ro')
%xlim([-15,5])
%ylim([-5,5])
xlabel('x (mm)','FontSize',16)
ylabel('y (mm)','FontSize',16)
figure;
hold on
contour(data.y,data.x,data.v,'b','ShowText','on')
contour(data.y,data.x,data.vT,'.-k','ShowText','on')
plot(0,0,'ro')
%xlim([-15,5])
%ylim([-5,5])
xlabel('x (mm)','FontSize',16)
ylabel('y (mm)','FontSize',16)
end

