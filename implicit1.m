clc;clear;


% 数据集的可变参数
alpha2 = 0;
alpha3 = 0.4;
alphaH = 0.6; %桁架的权重参数

t2 = 0.2;
t3 =-0.4;
lev=0.6; % 0.4~0.8 控制桁架体积比

X = linspace(-pi,pi,50);
Y = linspace(-pi,pi,50);
Z = linspace(-pi,pi,50);
[x,y,z] = meshgrid(X,Y,Z);


b1=1.3;
b2=0.5;
tPart=lev*(abs(x)+abs(y)+abs(z))/3+1;


x1=-pi;y1=-pi;z1=-pi;
x2=pi;y2=pi;z2=pi;

x0=(x1+x2)/2;y0=(y1+y2)/2;z0=(z1+z2)/2;LPart=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
dx2=x2-x0;dy2=y2-y0;dz2=z2-z0;
dx=x-x0;dy=y-y0;dz=z-z0;
Ld = sqrt(dx.^2+dy.^2+dz.^2);
cos_ThetaPart = sqrt(((dx*dx2+dy*dy2+dz*dz2)./(Ld*sqrt(dx2^2+dy2^2+dz2^2))).^2);
sin_ThetaPart = sqrt(1-cos_ThetaPart.^2);

FaiPartc1 =  (LPart/2)^2-(cos_ThetaPart.*Ld).^2;
FaiPartc2 =  (tPart/2).^2-(sin_ThetaPart.*Ld).^2;


FaiPartc_l1 = FaiPartc2;
%% 

% tPart=3;
x1=pi;y1=-pi;z1=-pi;
x2=-pi;y2=pi;z2=pi;

x0=(x1+x2)/2;y0=(y1+y2)/2;z0=(z1+z2)/2;LPart=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
dx2=x2-x0;dy2=y2-y0;dz2=z2-z0;
dx=x-x0;dy=y-y0;dz=z-z0;
Ld = sqrt(dx.^2+dy.^2+dz.^2);
cos_ThetaPart = sqrt(((dx*dx2+dy*dy2+dz*dz2)./(Ld*sqrt(dx2^2+dy2^2+dz2^2))).^2);
sin_ThetaPart = sqrt(1-cos_ThetaPart.^2);

FaiPartc1 =  (LPart/2)^2-(cos_ThetaPart.*Ld).^2;
FaiPartc2 =  (tPart/2).^2-(sin_ThetaPart.*Ld).^2;

FaiPartc_l2 = FaiPartc2;
%% 

% tPart=3;
x1=pi;y1=pi;z1=-pi;
x2=-pi;y2=-pi;z2=pi;

x0=(x1+x2)/2;y0=(y1+y2)/2;z0=(z1+z2)/2;LPart=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
dx2=x2-x0;dy2=y2-y0;dz2=z2-z0;
dx=x-x0;dy=y-y0;dz=z-z0;
Ld = sqrt(dx.^2+dy.^2+dz.^2);
cos_ThetaPart = sqrt(((dx*dx2+dy*dy2+dz*dz2)./(Ld*sqrt(dx2^2+dy2^2+dz2^2))).^2);
sin_ThetaPart = sqrt(1-cos_ThetaPart.^2);

FaiPartc1 =  (LPart/2)^2-(cos_ThetaPart.*Ld).^2;
FaiPartc2 =  (tPart/2).^2-(sin_ThetaPart.*Ld).^2;

FaiPartc_l3 = FaiPartc2;

%% 

% tPart=3;
x1=-pi;y1=pi;z1=-pi;
x2=pi;y2=-pi;z2=pi;

x0=(x1+x2)/2;y0=(y1+y2)/2;z0=(z1+z2)/2;LPart=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
dx2=x2-x0;dy2=y2-y0;dz2=z2-z0;
dx=x-x0;dy=y-y0;dz=z-z0;
Ld = sqrt(dx.^2+dy.^2+dz.^2);
cos_ThetaPart = sqrt(((dx*dx2+dy*dy2+dz*dz2)./(Ld*sqrt(dx2^2+dy2^2+dz2^2))).^2);
sin_ThetaPart = sqrt(1-cos_ThetaPart.^2);

FaiPartc1 =  (LPart/2)^2-(cos_ThetaPart.*Ld).^2;
FaiPartc2 =  (tPart/2).^2-(sin_ThetaPart.*Ld).^2;

FaiPartc_l4 = FaiPartc2;


FaiPartc_com = max(FaiPartc_l1,FaiPartc_l2);
FaiPartc_com = max(FaiPartc_l3,FaiPartc_com);
FaiPartc_com = max(FaiPartc_l4,FaiPartc_com)+lev;
% isosurface(x,y,z,FaiPartc_com,0)

figure(5001)
isosurface(x,y,z,FaiPartc_com,0)
isocaps(x,y,z,FaiPartc_com,0)




%%
alpha1 = 0;
% alpha2 = 0.8;
% alpha3 = 0.2;
t1 = 0;
% t2 = 0;
% t3 = 0;

% X = [-0.5:0.1:0.5]*2*pi;
% Y = [-0.5:0.1:0.5]*2*pi;
% Z = [-0.5:0.1:0.5]*2*pi;
% [x,y,z] = meshgrid(X,Y,Z);

F1 = cos(x) + cos(y) + cos(z)+ t1;
F2 = cos(x).*cos(y).*cos(z) - sin(x).*sin(y).*sin(z) + t2;
F3 = 8*cos(x).*cos(y).*cos(z)+cos(2*x).*cos(2*y).*cos(2*z)-(cos(2*x).*cos(2*y)+cos(2*y).*cos(2*z)+cos(2*z).*cos(2*x))+t3;
Fmerge = alpha1*4*F1+alpha2*4*F2+alpha3*F3;

figure(5002)
isosurface(x,y,z,Fmerge,0)



%% 
figure(5003)
isosurface(x,y,z,Fmerge+alphaH*FaiPartc_com,0)
isocaps(x,y,z,Fmerge+alphaH*FaiPartc_com,0)

xlabel('x')
ylabel('y')
zlabel('z')

% %% 转换为逻辑矩阵
% logic = zeros(50,50,50);
% 
% LevelSet = a1*Fmerge+a2*FaiPartc_com;
% logic(LevelSet>0) = 1;
% 
% CellStructrue(i).logic = logic;
