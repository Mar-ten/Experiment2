%% Mixed Mode Analysis
%filename={'Mixed_Mode-0010.csv','Mixed_Mode-0011.csv','Mixed_Mode-0012.csv','Mixed_Mode-0013.csv','Mixed_Mode-0014.csv','Mixed_Mode-0015.csv','Mixed_Mode-0016.csv','Mixed_Mode-0017.csv'};
filename={'Mode1-0010.csv','Mode1-0011.csv','Mode1-0012.csv','Mode1-0013.csv','Mode1-0014.csv','Mode1-0015.csv','Mode1-0016.csv','Mode1-0017.csv'};


nu=0.35; %Poisson's Ratio
E=3e3; % Elastic Modulus (MPa)
mu=E/(2*(1+nu)); % Shear Modulus from E and nu (MPa)
kappa=(3-nu)/(1+nu); % Plain Stress value for kappa
Load=[7.784,93.77,215.6,295,412,489.6,587,745,834,899,952,1010];
K = [0.0033,0.0392,0.0900,0.1232,0.1720,0.2044,0.2451,0.3111,0.3483,0.3754,0.3975,0.4218];
a=70; b=73;c=54;
K=K(end-length(filename)+1:end);
Load=Load(end-length(filename)+1:end);
Kexp1=zeros(size(K));
Kexp2=zeros(size(K));

for i=1:length(filename)
    raw=importdata(filename{i}); %import raw data
    data=[]; %create empty structure
    data.x=raw.data(:,2); %load data into structure
    data.y=raw.data(:,1);
    data.u=raw.data(:,4);
    data.v=raw.data(:,3);
    numY=sum(raw.data(:,5)==raw.data(1,5));  %calculate analysis area length and width
    numX=sum(raw.data(:,6)==raw.data(1,6));
    data.x=reshape(data.x,[numX,numY])/(1E3); %convert vectors to matricies and scale to meters
    data.y=reshape(data.y,[numX,numY])/(1E3);
    data.u=reshape(data.u,[numX,numY])/(1E3);
    data.v=reshape(data.v,[numX,numY])/(1E3);
    data.k_u=zeros(numX,numY); %create empty arrays for later
    data.k_u=zeros(numX,numY);
    data.uT=zeros(numX,numY);
    data.vT=zeros(numX,numY);
    F_1_u=zeros(numX,numY);
    F_1_v=zeros(numX,numY);
    F_2_u=zeros(numX,numY);
    F_2_v=zeros(numX,numY);
    K1=zeros(numX,numY);
    K2=zeros(numX,numY);
    for j=1:numX
        for k=1:numY
            r=sqrt(data.x(j,k)^2+data.y(j,k)^2);
            theta=atan2(data.y(j,k),data.x(j,k));
            theta(isnan(theta))=0;
            F_1_u=1/(8*mu*pi)*sqrt(2*pi*r)*((2*kappa-1)*cos(theta/2)-cos(3*theta/2));
            F_1_v=1/(8*mu*pi)*sqrt(2*pi*r)*((2*kappa+1)*sin(theta/2)-sin(3*theta/2));
            F_2_u=1/(8*mu*pi)*sqrt(2*pi*r)*((2*kappa+3)*sin(theta/2)+sin(3*theta/2));
            F_2_v=1/(8*mu*pi)*sqrt(2*pi*r)*((2*kappa-3)*cos(theta/2)+cos(3*theta/2));
            F=[F_1_u, F_1_v; F_2_u, F_2_v]';
            u=[data.u(j,k);data.v(j,k)];
            %             u=data.v(j,k);
            %             F=F_1_u;
            K_exp=F\u;
            K1(j,k)=K_exp(1);
            K2(j,k)=K_exp(2);
            
            
            %             data.k_u(j,k)=data.u(j,k)*(8*mu*pi)/(sqrt(2*pi*r)/((2*kappa-1)*cos(theta/2)-cos(3*theta/2))); %calculate experimental K_1
            %             data.k_v(j,k)=data.v(j,k)*(8*mu*pi)/(sqrt(2*pi*r)*((2*kappa+1)*sin(theta/2)-sin(3*theta/2)));
            %             data.vT(j,k)=K(i)/(8*mu*pi)*sqrt(2*pi*r)*((2*kappa-1)*cos(theta/2)-cos(3*theta/2)); %calculate theoretical displacements
            %             data.uT(j,k)=K(i)/(8*mu*pi)*sqrt(2*pi*r)*((2*kappa+1)*sin(theta/2)-sin(3*theta/2));
        end
    end
    %
    %     data.k_u(abs(data.k_u)>=2)=0; % gid rid of large values where denominator approaches 0
    %     data.k_v(abs(data.k_v)>=2)=0; % gid rid of large values where denominator approaches 0
    %     Kexp1(i)=abs(data.k_v(a,b)); % Record experimental K_1 value at point 1
    %     Kexp2(i)=abs(data.k_v(c,b)); % Record experimental K_1 value at point 2
    
    K1=abs(K1); K1(K1>=1)=1;
    K2=abs(K2); K2(K2>=1)=1;
    data.k1=K1;
    data.k2=K2;
%     assignin('base',['raw',num2str(i)],raw); %create raw data structure to access outside the loop
    assignin('base',['data',num2str(i)],data); %%create processed data structure to access outside the loop
end

% % Rsq1=round(corr(K',Kexp1','type','Pearson')^2,2);
% % Rsq2=round(corr(K',Kexp2','type','Pearson')^2,2);
% Rsq1=corr(K',Kexp1','type','Pearson')^2;
% Rsq2=corr(K',Kexp2','type','Pearson')^2;
% figure; % Plot Experimental vs. Theoretical K_1 Value
% hold on
% plot(Load,K)
% plot(Load,Kexp1,'kd')
% plot(Load,Kexp2,'ro')
% xlabel('$\textrm{Load} (N)$','FontSize',16,'Interpreter','latex')
% ylabel('$K_1 (MPa \sqrt{m})$','FontSize',16,'Interpreter','latex')
% legend({'Theoretical K_1',['Experimental K_1 Pt1, ','R^2=',num2str(Rsq1)],['Experimental K_1 Pt2, ','R^2=',num2str(Rsq2)]},'Location','NorthWest','FontSize',11)

%%

figure; contourf(data8.x*1E3,data8.y*1E3,data8.k1); colorbar
xlabel('x (mm)','FontSize',16)
ylabel('y (mm)','FontSize',16)
title('K_1','FontSize',20)
colorbar
figure; contourf(data8.x*1E3,data8.y*1E3,data8.k2); colorbar
xlabel('x (mm)','FontSize',16)
ylabel('y (mm)','FontSize',16)
title('K_2','FontSize',20)
colorbar

ratio=data8.k1./data8.k2;
ratio(ratio>=10)=10;
figure; contourf(data8.x,data8.y,ratio); colorbar
