
function [Xcm,Ycm]=get_CM(img,thr)
%thr is the threshold below which the program consider it zero
I=rescale(img);
I(I<thr)=0;
X=1:size(I,2);
Y=1:size(I,1);
[Xg, Yg]= meshgrid(X,Y);
TotalMass=sum(I,'all');
Xcm=sum(I.*Xg,'all')/(TotalMass);%*size(I,2)); not normalized to the sensor size
Ycm=sum(I.*Yg,'all')/(TotalMass);%*size(I,1)); %not normalized to the sensor size

end
