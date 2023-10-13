clc;
clear;

rng default % For reproducibility
X = lhsdesign(1000,7);


for i = 1:1000
    Y(i,1) = X(i,1)/(X(i,1)+X(i,2)+X(i,3));
    Y(i,2) = X(i,2)/(X(i,1)+X(i,2)+X(i,3));
    Y(i,3) = X(i,3)/(X(i,1)+X(i,2)+X(i,3));
end

T1 = X(:,4)*0.8-0.4;
T2 = X(:,5)*0.8-0.4;
T3 = X(:,6)*1.3-0.65;
U = X(:,7)*20+1;

LaStrSet = [Y T1 T2 T3 U];
% 数据集的可变参数
tPart = 1.5;
alpha1 = 0;
alpha2 = 0.1;
alpha3 = 0.1;
alphaH = 0.8;
t2 = 0.5;
t3 = 0.5;
t4 = 0.8; % 范围可以取-0.3~1
u = 5; % 范围可以取1~21之间

alpha7 = 0.5;
alpha8 = 0.5;


FaiPartc_com1 = 0;
FaiPartc_com2 = 0;
X = linspace(-pi,pi,50);
Y = linspace(-pi,pi,50);
Z = linspace(-pi,pi,50);
[x,y,z] = meshgrid(X,Y,Z);

%% 全立方，缺
endpoint1 = [% 八斜边
    0 0 pi 0 -pi 0
    0 0 pi pi 0 0
    0 0 pi 0 pi 0
    0 0 pi -pi 0 0
    0 0 -pi 0 -pi 0
    0 0 -pi pi 0 0
    0 0 -pi 0 pi 0
    0 0 -pi -pi 0 0
    % 中心三立柱
    pi 0 0 -pi 0 0
    0 pi 0 0 -pi 0
    0 0 pi 0 0 -pi
    % x轴垂面
    pi -pi pi pi pi -pi
    pi pi pi pi -pi -pi
    -pi -pi pi -pi pi -pi
    -pi pi pi -pi -pi -pi
    % y垂面
    pi pi pi -pi pi -pi
    -pi pi pi pi pi -pi
    pi -pi pi -pi -pi -pi
    -pi -pi pi pi -pi -pi
    % z垂面
    pi -pi pi -pi pi pi
    pi pi pi -pi -pi pi
    pi -pi -pi -pi pi -pi
    pi pi -pi -pi -pi -pi
    % 周面立柱
    pi -pi pi pi -pi -pi
    pi pi pi pi pi -pi
    -pi pi pi -pi pi -pi
    -pi -pi pi -pi -pi -pi
    % 顶面立柱
    pi -pi pi pi pi pi
    pi pi pi -pi pi pi
    -pi pi pi -pi -pi pi
    -pi -pi pi pi -pi pi
    % 底面立柱
    pi -pi -pi pi pi -pi
    pi pi -pi -pi pi -pi
    -pi pi -pi -pi -pi -pi
    -pi -pi -pi pi -pi -pi
    % 补充四斜边
        pi 0 0 0 pi 0
    0 pi 0 -pi 0 0 
    -pi 0 0 0 -pi 0
    0 -pi 0 pi 0 0
    ];
%% 八角
endpoint1 = [
    -pi -pi -pi pi pi pi
    pi -pi -pi -pi pi pi
    pi pi -pi -pi -pi pi
    -pi pi -pi pi -pi pi
    ];

%% 六面交叉
endpoint1 = [
    pi -pi pi -pi pi pi
    -pi -pi pi pi pi pi
    -pi -pi -pi pi pi -pi
    pi -pi -pi -pi pi -pi
    -pi -pi pi pi -pi -pi
    pi -pi pi -pi -pi -pi
    -pi pi pi pi pi -pi
    pi pi pi -pi pi -pi
    pi -pi pi pi pi -pi
    pi -pi -pi pi pi pi
    -pi pi pi -pi -pi -pi
    -pi -pi pi -pi pi -pi
    ];

%% cube

endpoint1 = [
    pi -pi pi pi pi pi
    pi pi pi -pi pi pi
    -pi pi pi -pi -pi pi
    pi -pi pi -pi -pi pi
    pi -pi -pi pi pi -pi 
    pi pi -pi -pi pi -pi
     -pi -pi -pi -pi pi -pi
     -pi -pi -pi pi -pi -pi
     pi -pi pi pi -pi -pi
     pi pi pi pi pi -pi 
     -pi pi pi -pi pi -pi
      -pi -pi pi -pi -pi -pi
    ];

%% 六面交叉++
endpoint1 = [
    pi -pi pi -pi pi pi
    -pi -pi pi pi pi pi
    -pi -pi -pi pi pi -pi
    pi -pi -pi -pi pi -pi
    -pi -pi pi pi -pi -pi
    pi -pi pi -pi -pi -pi
    -pi pi pi pi pi -pi
    pi pi pi -pi pi -pi
    pi -pi pi pi pi -pi
    pi -pi -pi pi pi pi
    -pi pi pi -pi -pi -pi
    -pi -pi pi -pi pi -pi
    % 八斜边
    0 0 pi 0 -pi 0
    0 0 pi pi 0 0
    0 0 pi 0 pi 0
    0 0 pi -pi 0 0
    0 0 -pi 0 -pi 0
    0 0 -pi pi 0 0
    0 0 -pi 0 pi 0
    0 0 -pi -pi 0 0
    pi 0 0 0 pi 0
    0 pi 0 -pi 0 0 
    -pi 0 0 0 -pi 0
    0 -pi 0 pi 0 0
    ];

%% 三立柱
endpoint1 = [
    0 0 pi 0 0 -pi
    0 -pi 0 0 pi 0
    pi 0 0 -pi 0 0
    ];

%% 三立柱
endpoint1 = [
    0 0 pi 0 0 -pi
    0 -pi 0 0 pi 0
    pi 0 0 -pi 0 0
    ];

%% 三立柱+八角
endpoint1 = [
    0 0 pi 0 0 -pi
    0 -pi 0 0 pi 0
    pi 0 0 -pi 0 0
    -pi -pi -pi pi pi pi
    pi -pi -pi -pi pi pi
    pi pi -pi -pi -pi pi
    -pi pi -pi pi -pi pi
    ];
%% 六面交叉+ 中心斜边
endpoint1 = [
    pi -pi pi -pi pi pi
    -pi -pi pi pi pi pi
    -pi -pi -pi pi pi -pi
    pi -pi -pi -pi pi -pi
    -pi -pi pi pi -pi -pi
    pi -pi pi -pi -pi -pi
    -pi pi pi pi pi -pi
    pi pi pi -pi pi -pi
    pi -pi pi pi pi -pi
    pi -pi -pi pi pi pi
    -pi pi pi -pi -pi -pi
    -pi -pi pi -pi pi -pi

    ];
%% 六面交叉+ cube 立柱
endpoint1 = [
    pi -pi pi -pi pi pi
    -pi -pi pi pi pi pi
    -pi -pi -pi pi pi -pi
    pi -pi -pi -pi pi -pi
    -pi -pi pi pi -pi -pi
    pi -pi pi -pi -pi -pi
    -pi pi pi pi pi -pi
    pi pi pi -pi pi -pi
    pi -pi pi pi pi -pi
    pi -pi -pi pi pi pi
    -pi pi pi -pi -pi -pi
    -pi -pi pi -pi pi -pi
        pi -pi pi pi pi pi
    pi pi pi -pi pi pi
    -pi pi pi -pi -pi pi
    pi -pi pi -pi -pi pi
    pi -pi -pi pi pi -pi 
    pi pi -pi -pi pi -pi
     -pi -pi -pi -pi pi -pi
     -pi -pi -pi pi -pi -pi
     pi -pi pi pi -pi -pi
     pi pi pi pi pi -pi 
     -pi pi pi -pi pi -pi
      -pi -pi pi -pi -pi -pi
           pi -pi pi pi pi pi
    pi pi pi -pi pi pi
    -pi pi pi -pi -pi pi
    pi -pi pi -pi -pi pi
    pi -pi -pi pi pi -pi 
    pi pi -pi -pi pi -pi
     -pi -pi -pi -pi pi -pi
     -pi -pi -pi pi -pi -pi
     pi -pi pi pi -pi -pi
     pi pi pi pi pi -pi 
     -pi pi pi -pi pi -pi
      -pi -pi pi -pi -pi -pi
    ];
%% 六面交叉
endpoint1 = pi*0.1*[
    -pi -pi pi pi -pi -pi
    pi -pi pi -pi -pi -pi
    -pi pi pi pi pi -pi
    pi pi pi -pi pi -pi
    pi -pi pi pi pi -pi
    pi -pi -pi pi pi pi
    -pi pi pi -pi -pi -pi
    -pi -pi pi -pi pi -pi
    ]/pi;
% %%
%     pd = makedist('Normal');
% endpoint1 = random(pd,[4,6]);
% endpoint2 = endpoint1;
% endpoint2(:,1) = -endpoint1(:,1);
% endpoint2(:,4) = -endpoint1(:,4);
% endpoint1 = [endpoint1
%     endpoint2];
% 
% endpoint2 = endpoint1;
% endpoint2(:,2) = -endpoint1(:,2);
% endpoint2(:,5) = -endpoint1(:,5);
% endpoint1 = [endpoint1
%     endpoint2];
% endpoint2 = endpoint1;
% endpoint2(:,3) = -endpoint1(:,3);
% endpoint2(:,6) = -endpoint1(:,6);
% endpoint1 = [endpoint1
%     endpoint2];




for j = 1:size(LaStrSet,1)

    FaiPartc_com1 = 0;

 

    for i =1:size(endpoint1,1)
        FaiPartc_com1 = StrutGen(endpoint1(i,:),1,FaiPartc_com1);
    end

    % for i =1:size(endpoint2,1)
    %     FaiPartc_com2 = StrutGen(endpoint2(i,:),1.2,FaiPartc_com2);
    % end

figure(5001)
clf(5001)
isosurface(x,y,z,FaiPartc_com1,0)
isocaps(x,y,z,FaiPartc_com1,0)

% figure(5004)
% clf(5004)
% isosurface(x,y,z,FaiPartc_com2,0)
% isocaps(x,y,z,FaiPartc_com2,0)

%
% alpha1 = 0;
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
Fmerge = alpha1*4*F1+alpha2*F2+alpha3*F3;

figure(5002)
clf(5002)
isosurface(x,y,z,Fmerge,0)
isocaps(x,y,z,Fmerge,0)

figure(5006)
clf(5006)
isosurface(x,y,z,Fmerge+alphaH*4*FaiPartc_com1,0)
isocaps(x,y,z,Fmerge+alphaH*4*FaiPartc_com1,0)

% figure(5005)
% isosurface(x,y,z,alpha7*FaiPartc_com1+alpha8*FaiPartc_com2,0)
% isocaps(x,y,z,alpha7*FaiPartc_com1+alpha8*FaiPartc_com2,0)

%% 导出逻辑矩阵

logic = zeros(50,50,50);

LevelSet = Fmerge+alphaH*FaiPartc_com1;
logic(LevelSet>0) = 1;

CellStructrue(j).logic = logic;

FaiPartc_com = 0;

xlabel('x')
ylabel('y')
zlabel('z')

end


%% 支柱水平集矩阵生成函数

function FaiPartc_com = StrutGen(endpoint,tPart,FaiPartc_com)

X = linspace(-pi,pi,50);
Y = linspace(-pi,pi,50);
Z = linspace(-pi,pi,50);
[x,y,z] = meshgrid(X,Y,Z);

x1=endpoint(1);y1=endpoint(2);z1=endpoint(3);
x2=endpoint(4);y2=endpoint(5);z2=endpoint(6);

x0=(x1+x2)/2;y0=(y1+y2)/2;z0=(z1+z2)/2;LPart=sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
dx2=x2-x0;dy2=y2-y0;dz2=z2-z0;
dx=x-x0;dy=y-y0;dz=z-z0;
Ld = sqrt(dx.^2+dy.^2+dz.^2);
cos_ThetaPart = sqrt(((dx*dx2+dy*dy2+dz*dz2)./(Ld*sqrt(dx2^2+dy2^2+dz2^2))).^2);
sin_ThetaPart = sqrt(1-cos_ThetaPart.^2);

FaiPartc1 =  (LPart/2)^2-(cos_ThetaPart.*Ld).^2;
FaiPartc2 =  (tPart/2).^2-(sin_ThetaPart.*Ld).^2;


FaiPartc_l1 = FaiPartc2;

if FaiPartc_com == 0
    FaiPartc_com = FaiPartc_l1;
else
    FaiPartc_com = max(FaiPartc_l1,FaiPartc_com);
end


end
