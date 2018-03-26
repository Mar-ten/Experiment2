data=importdata('35 fracture-0017_0.csv');
%close all
x=data.data(:,1);
y=data.data(:,2);
numY=sum(data.data(:,5)==data.data(1,5));
numX=sum(data.data(:,6)==data.data(1,6));


nu=0.35;
E=3e3;
mu=E/(2*(1+nu));
kappa=(3-nu)/(1+nu);

u=data.data(:,3);
v=data.data(:,4);
x2=reshape(x,[numX,numY]);
y2=reshape(y,[numX,numY]);
u2=reshape(u,[numX,numY]);
v2=reshape(v,[numX,numY]);

k_u=zeros(numX,numY);
k_v=zeros(numX,numY);
r_calc=zeros(numX,numY);
angle=zeros(numX,numY);

x3=x2/(1E3); %convert to meters
y3=y2/(1E3); %convert to meters
u3=u2/(1E3); %convert to meters
v3=v2/(1E3); %convert to meters

for i=1:numX
    for j=1:numY
        r=sqrt(x3(i,j)^2+y3(i,j)^2);
        theta=atan2(x3(i,j),y3(i,j));
        theta(isnan(theta))=0;
        angle(i,j)=theta;
        r_calc(i,j)=r;
        k_u(i,j)=u3(i,j)*(8*mu*pi)/(sqrt(2*pi*r)/((2*kappa-1)*cos(theta/2)-cos(3*theta/2)));
        k_v(i,j)=v3(i,j)*(8*mu*pi)/(sqrt(2*pi*r)*((2*kappa+1)*sin(theta/2)-sin(3*theta/2)));
        
    end
end

figure; contourf(y3,x3,u3)
figure; contourf(y3,x3,v3)
figure;quiver(y3,x3,v3,u3)

k_u(abs(k_u)>=2)=0;
figure; contourf(y3,x3,abs(k_u))

k_v(abs(k_v)>=2)=0;
figure; contourf(y3,x3,abs(k_v))