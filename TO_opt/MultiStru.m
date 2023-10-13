clc;
clear;

rng default % For reproducibility
X = lhsdesign(1000,5);
X = X*3;

X(X<0.5) = 0;

% 数据集的可变参数
tpart1 = X(:,1);
tpart2 = X(:,2);
tpart3 = X(:,3);
tpart4 = X(:,4);
tpart5 = X(:,5);



FaiPartc_com1 = 0;
FaiPartc_com2 = 0;
X = linspace(-pi,pi,50);
Y = linspace(-pi,pi,50);
Z = linspace(-pi,pi,50);
[x,y,z] = meshgrid(X,Y,Z);

% %% 全立方，缺
% endpoint1 = [% 八斜边
%     0 0 pi 0 -pi 0
%     0 0 pi pi 0 0
%     0 0 pi 0 pi 0
%     0 0 pi -pi 0 0
%     0 0 -pi 0 -pi 0
%     0 0 -pi pi 0 0
%     0 0 -pi 0 pi 0
%     0 0 -pi -pi 0 0
%     % 中心三立柱
%     pi 0 0 -pi 0 0
%     0 pi 0 0 -pi 0
%     0 0 pi 0 0 -pi
%     % x轴垂面
%     pi -pi pi pi pi -pi
%     pi pi pi pi -pi -pi
%     -pi -pi pi -pi pi -pi
%     -pi pi pi -pi -pi -pi
%     % y垂面
%     pi pi pi -pi pi -pi
%     -pi pi pi pi pi -pi
%     pi -pi pi -pi -pi -pi
%     -pi -pi pi pi -pi -pi
%     % z垂面
%     pi -pi pi -pi pi pi
%     pi pi pi -pi -pi pi
%     pi -pi -pi -pi pi -pi
%     pi pi -pi -pi -pi -pi
%     % 周面立柱
%     pi -pi pi pi -pi -pi
%     pi pi pi pi pi -pi
%     -pi pi pi -pi pi -pi
%     -pi -pi pi -pi -pi -pi
%     % 顶面立柱
%     pi -pi pi pi pi pi
%     pi pi pi -pi pi pi
%     -pi pi pi -pi -pi pi
%     -pi -pi pi pi -pi pi
%     % 底面立柱
%     pi -pi -pi pi pi -pi
%     pi pi -pi -pi pi -pi
%     -pi pi -pi -pi -pi -pi
%     -pi -pi -pi pi -pi -pi
%     % 补充四斜边
%         pi 0 0 0 pi 0
%     0 pi 0 -pi 0 0 
%     -pi 0 0 0 -pi 0
%     0 -pi 0 pi 0 0
%     ];
% %% 八角
% endpoint1 = [
%     -pi -pi -pi pi pi pi
%     pi -pi -pi -pi pi pi
%     pi pi -pi -pi -pi pi
%     -pi pi -pi pi -pi pi
%     ];
% 
% %% 六面交叉
% endpoint1 = [
%     pi -pi pi -pi pi pi
%     -pi -pi pi pi pi pi
%     -pi -pi -pi pi pi -pi
%     pi -pi -pi -pi pi -pi
%     -pi -pi pi pi -pi -pi
%     pi -pi pi -pi -pi -pi
%     -pi pi pi pi pi -pi
%     pi pi pi -pi pi -pi
%     pi -pi pi pi pi -pi
%     pi -pi -pi pi pi pi
%     -pi pi pi -pi -pi -pi
%     -pi -pi pi -pi pi -pi
%     ];
% 
% %% cube
% 
% endpoint1 = [
%     pi -pi pi pi pi pi
%     pi pi pi -pi pi pi
%     -pi pi pi -pi -pi pi
%     pi -pi pi -pi -pi pi
%     pi -pi -pi pi pi -pi 
%     pi pi -pi -pi pi -pi
%      -pi -pi -pi -pi pi -pi
%      -pi -pi -pi pi -pi -pi
%      pi -pi pi pi -pi -pi
%      pi pi pi pi pi -pi 
%      -pi pi pi -pi pi -pi
%       -pi -pi pi -pi -pi -pi
%     ];
% 
% %% 六面交叉++
% endpoint1 = [
%     pi -pi pi -pi pi pi
%     -pi -pi pi pi pi pi
%     -pi -pi -pi pi pi -pi
%     pi -pi -pi -pi pi -pi
%     -pi -pi pi pi -pi -pi
%     pi -pi pi -pi -pi -pi
%     -pi pi pi pi pi -pi
%     pi pi pi -pi pi -pi
%     pi -pi pi pi pi -pi
%     pi -pi -pi pi pi pi
%     -pi pi pi -pi -pi -pi
%     -pi -pi pi -pi pi -pi
%     % 八斜边
%     0 0 pi 0 -pi 0
%     0 0 pi pi 0 0
%     0 0 pi 0 pi 0
%     0 0 pi -pi 0 0
%     0 0 -pi 0 -pi 0
%     0 0 -pi pi 0 0
%     0 0 -pi 0 pi 0
%     0 0 -pi -pi 0 0
%     pi 0 0 0 pi 0
%     0 pi 0 -pi 0 0 
%     -pi 0 0 0 -pi 0
%     0 -pi 0 pi 0 0
%     ];
% 
% %% 三立柱
% endpoint1 = [
%     0 0 pi 0 0 -pi
%     0 -pi 0 0 pi 0
%     pi 0 0 -pi 0 0
%     ];
% 
% %% 三立柱
% endpoint1 = [
%     0 0 pi 0 0 -pi
%     0 -pi 0 0 pi 0
%     pi 0 0 -pi 0 0
%     ];
% 
% %% 三立柱+八角
% endpoint1 = [
%     0 0 pi 0 0 -pi
%     0 -pi 0 0 pi 0
%     pi 0 0 -pi 0 0
%     -pi -pi -pi pi pi pi
%     pi -pi -pi -pi pi pi
%     pi pi -pi -pi -pi pi
%     -pi pi -pi pi -pi pi
%     ];
% %% 六面交叉+ 中心斜边
% endpoint1 = [
%     pi -pi pi -pi pi pi
%     -pi -pi pi pi pi pi
%     -pi -pi -pi pi pi -pi
%     pi -pi -pi -pi pi -pi
%     -pi -pi pi pi -pi -pi
%     pi -pi pi -pi -pi -pi
%     -pi pi pi pi pi -pi
%     pi pi pi -pi pi -pi
%     pi -pi pi pi pi -pi
%     pi -pi -pi pi pi pi
%     -pi pi pi -pi -pi -pi
%     -pi -pi pi -pi pi -pi
% 
%     ];
% %% 六面交叉+ cube 立柱
% endpoint2 = [
%     pi -pi pi -pi pi pi
%     -pi -pi pi pi pi pi
%     -pi -pi -pi pi pi -pi
%     pi -pi -pi -pi pi -pi
%     -pi -pi pi pi -pi -pi
%     pi -pi pi -pi -pi -pi
%     -pi pi pi pi pi -pi
%     pi pi pi -pi pi -pi
%     pi -pi pi pi pi -pi
%     pi -pi -pi pi pi pi
%     -pi pi pi -pi -pi -pi
%     -pi -pi pi -pi pi -pi
%         pi -pi pi pi pi pi
%     pi pi pi -pi pi pi
%     -pi pi pi -pi -pi pi
%     pi -pi pi -pi -pi pi
%     pi -pi -pi pi pi -pi 
%     pi pi -pi -pi pi -pi
%      -pi -pi -pi -pi pi -pi
%      -pi -pi -pi pi -pi -pi
%      pi -pi pi pi -pi -pi
%      pi pi pi pi pi -pi 
%      -pi pi pi -pi pi -pi
%       -pi -pi pi -pi -pi -pi
%            pi -pi pi pi pi pi
%     pi pi pi -pi pi pi
%     -pi pi pi -pi -pi pi
%     pi -pi pi -pi -pi pi
%     pi -pi -pi pi pi -pi 
%     pi pi -pi -pi pi -pi
%      -pi -pi -pi -pi pi -pi
%      -pi -pi -pi pi -pi -pi
%      pi -pi pi pi -pi -pi
%      pi pi pi pi pi -pi 
%      -pi pi pi -pi pi -pi
%       -pi -pi pi -pi -pi -pi
%     ];
% %% 六面交叉
% endpoint1 = pi*0.1*[
%     -pi -pi pi pi -pi -pi
%     pi -pi pi -pi -pi -pi
%     -pi pi pi pi pi -pi
%     pi pi pi -pi pi -pi
%     pi -pi pi pi pi -pi
%     pi -pi -pi pi pi pi
%     -pi pi pi -pi -pi -pi
%     -pi -pi pi -pi pi -pi
%     ]/pi;

%% 三立柱
endpoint1 = [
    0 0 pi 0 0 -pi
    0 -pi 0 0 pi 0
    pi 0 0 -pi 0 0
    ];
%% 中心斜柱
endpoint2 = [
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
%% 八角
endpoint3 = [
    -pi -pi -pi pi pi pi
    pi -pi -pi -pi pi pi
    pi pi -pi -pi -pi pi
    -pi pi -pi pi -pi pi
    ];
%% 六面交叉
endpoint4 = [
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
endpoint5 = [
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




for j = 1:size(tpart1,1)

    com1 = 0;
    com2 = 0;
    com3 = 0;
    com4 = 0;
    com5 = 0;
 

    for i =1:size(endpoint1,1)
        com1 = StrutGen(endpoint1(i,:),tpart1(j),com1);
    end

    for i =1:size(endpoint2,1)
        com2 = StrutGen(endpoint2(i,:),tpart2(j),com2);
    end

    for i =1:size(endpoint3,1)
        com3 = StrutGen(endpoint3(i,:),tpart3(j),com3);
    end
  
    for i =1:size(endpoint4,1)
        com4 = StrutGen(endpoint4(i,:),tpart4(j),com4);
    end
    
    for i =1:size(endpoint5,1)
        com5 = StrutGen(endpoint5(i,:),tpart5(j),com5);
    end


% figure(5001)
% clf(5001)
% isosurface(x,y,z,com1,0)
% isocaps(x,y,z,com1,0)
% 
% figure(5004)
% clf(5004)
% isosurface(x,y,z,com2,0)
% isocaps(x,y,z,com2,0)

ComCom = zeros(50,50,50,5);
ComCom(:,:,:,1) = com1;
ComCom(:,:,:,2) = com2;
ComCom(:,:,:,3) = com3;
ComCom(:,:,:,4) = com4;
ComCom(:,:,:,5) = com5;


a = max(ComCom,[],4);

figure(5006)
clf(5006)
isosurface(x,y,z,a,0)
isocaps(x,y,z,a,0)

% figure(5005)
% isosurface(x,y,z,alpha7*FaiPartc_com1+alpha8*FaiPartc_com2,0)
% isocaps(x,y,z,alpha7*FaiPartc_com1+alpha8*FaiPartc_com2,0)

%% 导出逻辑矩阵

logic = zeros(50,50,50);

LevelSet = a;
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
