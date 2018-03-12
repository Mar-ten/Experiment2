function [ U,V ] = DIC2IMG( file )
%DIC2IMG takes DIC data written to a CSV file and converts it to a matlab
%array for processing

% example: [ U,V ] = DIC2IMG( 'painted_mode_1_data.csv )

Data=importdata(file);
xmin=min(Data.data(:,3));
xmax=max(Data.data(:,3));
ymin=min(Data.data(:,4));
ymax=max(Data.data(:,4));

L=xmax-xmin;
H=ymax-ymin;

U=zeros(L,H);
V=zeros(L,H);

for i=1:length(Data.data)
    U(Data.data(i,3)-xmin+1,Data.data(i,4)-ymin+1)=Data.data(i,5);
    V(Data.data(i,3)-xmin+1,Data.data(i,4)-ymin+1)=Data.data(i,6);
end

U=flip(U'); %Changes Orientation of plots
V=flip(V');


figure;
imagesc(U); colorbar
figure;
imagesc(V); colorbar
end

