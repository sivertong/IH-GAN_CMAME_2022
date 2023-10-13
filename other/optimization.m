function optimization
% 定义材料参数
E0 = 210e9;           % 杨氏模量
nu = 0.3;             % 泊松比
% 加载像素模型并修正
load('x.mat');
load('x_v.mat');
x(x_v==3|x_v==2)=1;
nelx=size(x,2); % x方向单元数量
nely=size(x,1); % y方向单元数量
nelz=size(x,3); % z方向单元数量
volfrac=0.5;    % 预设体分比
nele = nelx*nely*nelz; % 单元总数
% 定义设计域
DesignD=find(x);
NumSolid=length(DesignD);     %初始模型中的实体单元数量
Tpercent=ceil(NumSolid*0.02); %每个循环去除单元数量
% 优化循环参数
maxloop = 200;    % 最大迭代次数
tolx = 0.01;      % 迭代终止条件
% 识别非设计单元，加载单元和固定单元
FixEleInd=zeros(size(x));    % 固定单元索引
passiveEleInd=zeros(size(x));% 非设计单元索引
LoadEleInd=zeros(size(x));   % 加载单元索引
for ix=1:nelx
    for iy=1:nely
        for iz=1:nelz
            if (ix>8&&ix<21&&iy>95&&iy<106&&iz>4&&iz<11)&&x(iy,ix,iz)==0
                LoadEleInd(iy,ix,iz)=1;
            end
            if (ix>107&&ix<116&&iy>109&&iy<123&&iz>4&&iz<11)&&x(iy,ix,iz)==0
                FixEleInd(iy,ix,iz)=2;
            end
            if iy<10&&x(iy,ix,iz)==1
                FixEleInd(iy,ix,iz)=1;
            end
            if (((ix-13.5)^2+(iy-100.5)^2-7.5^2<0)||(ix>106&&iy>108))||(iy<10&&x(iy,ix,iz)==1)
                passiveEleInd(iy,ix,iz)=1;
            end
        end
    end
end
% 过滤预处理
rmin=2;
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
% 定义结点编号和自由度
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
edofMat_new = edofMat(x(:)==1|passiveEleInd(:)==1,:);  %单元自由度编号矩阵，每一行代表一个单元的自由度编号
ndof_new = unique(edofMat_new(:));
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% 识别固定单元对应的自由度
FixedDofIndx1=edofMat(FixEleInd==1,[1:3:24,3:3:24]); %约束x方向和z方向位移
fixeddof1 = unique(FixedDofIndx1(:));
FixedDofIndx2=edofMat(FixEleInd==2,:);
fixeddof2 = unique(FixedDofIndx2(:));
FixedDofIndx3=edofMat(LoadEleInd==1,3:3:24);%约束z方向位移
fixeddof3 = unique(FixedDofIndx3(:));
fixeddof11=union(fixeddof1,fixeddof2);
fixeddof=union(fixeddof11,fixeddof3); %固定的自由度编号
% 识别加载单元对应的自由度
LoadDofIndx1=edofMat(LoadEleInd==1,2:3:24); %y方向施加载荷
loaddof1 = unique(LoadDofIndx1(:));
LoadDofIndx2=edofMat(LoadEleInd==1,1:3:24); %x方向施加载荷
loaddof2 = unique(LoadDofIndx2(:));
% 定义载荷矩阵和位移向量
F1 = sparse(loaddof1,1,8.737/length(loaddof1),ndof,1);
F2 = sparse(loaddof2,1,8.197/length(loaddof1),ndof,1);
F=[F1,F2];
U1 = zeros(ndof,1);
U2 = zeros(ndof,1);
freedofs = setdiff(ndof_new,fixeddof);
KE = lk_H8(nu); %单元刚度矩阵
% 初始化优化循环
loop = 0;
change = 1;
Record=[];
% 开始迭代
while change > tolx && loop < maxloop
    loop = loop+1;
    xpassive=x+9*passiveEleInd;
    sK = reshape(KE(:)*xpassive(:)',24*24*nele,1);
    K = sparse(iK,jK,sK);
    %% DCG 预处理
    lx = repmat((0:nelx)',1,nely+1,nelz+1);
    ly = repmat((0:nely),nelx+1,1,nelz+1);
    lz = repmat(reshape(0:nelz,1,1,nelz+1),nely+1,nelx+1,1);
    nel = (nelx+1)*(nely+1)*(nelz+1);
    W1 = repmat(diag(ones(3,1)),nel,1);
    W21 = [0,0,0;0,0,-1;0,1,0]*[lx(:)';ly(:)';lz(:)'];
    W22 = [0,0,1;0,0,0;-1,0,0]*[lx(:)';ly(:)';lz(:)'];
    W23 = [0,-1,0;1,0,0;0,0,0]*[lx(:)';ly(:)';lz(:)'];
    num_mat = reshape(1:nel,(nelx+1),(nely+1),(nelz+1));
    group_x_num = 1;
    group_y_num = 1;
    group_z_num = 1;
    group_x_begin = 1:floor(nelx/group_x_num):nelx;  group_x_end = floor(nelx/group_x_num):floor(nelx/group_x_num):nelx;  group_x_end(group_x_num)=nelx+1;
    group_y_begin = 1:floor(nely/group_y_num):nely;  group_y_end = floor(nely/group_y_num):floor(nely/group_y_num):nely;  group_y_end(group_y_num)=nely+1;
    group_z_begin = 1:floor(nelz/group_z_num):nelz;  group_z_end = floor(nelz/group_z_num):floor(nelz/group_z_num):nelz;  group_z_end(group_z_num)=nelz+1;
    group_x=cell(group_x_num,1); for i=1:group_x_num; group_x{i}=group_x_begin(i):group_x_end(i); end
    group_y=cell(group_y_num,1); for i=1:group_y_num; group_y{i}=group_y_begin(i):group_y_end(i); end
    group_z=cell(group_z_num,1); for i=1:group_z_num; group_z{i}=group_z_begin(i):group_z_end(i); end
    group_mat = cell(numel(group_x)*numel(group_y),1); mean_mat=group_mat; group_mat_dof=group_mat;
    mean1_mat=mean_mat; mean2_mat=mean_mat; mean3_mat=mean_mat;
    for k=1:numel(group_z)
        for i=1:numel(group_x)
            for j=1:numel(group_y)
                ni = j + (i-1)*numel(group_y) + (k-1)*numel(group_x)*numel(group_y);
                group_mat{ni} = num_mat(group_x{i},group_y{j},group_z{k});
                group_mat_dof{ni} = [3*group_mat{ni}(:)-2,3*group_mat{ni}(:)-1,3*group_mat{ni}(:)]';
                mean1_mat{ni} = mean(W21(:,group_mat{ni}(:)),2);
                mean2_mat{ni} = mean(W22(:,group_mat{ni}(:)),2);
                mean3_mat{ni} = mean(W23(:,group_mat{ni}(:)),2);
                W21(:,group_mat{ni}(:)) = W21(:,group_mat{ni}(:)) - mean1_mat{ni};
                W22(:,group_mat{ni}(:)) = W22(:,group_mat{ni}(:)) - mean2_mat{ni};
                W23(:,group_mat{ni}(:)) = W23(:,group_mat{ni}(:)) - mean3_mat{ni};
            end
        end
    end
    W = [W1,W21(:),W22(:),W23(:)];
    W_diag = [];
    for i=1:length(group_mat_dof)
        W_diag = blkdiag( W_diag , W(group_mat_dof{i}(:),:) );
    end
    %% 免组装
    a1=edofMat;
    a1(ismember(edofMat(:),freedofs)==0)=nan;
    a1_logcal = isfinite(a1);
    WKW_free=0;
    for i=1:length(a1)
        logcal_element = a1_logcal(i,:);
        if all(logcal_element==0); continue; end
        num_a1 = a1(i,a1_logcal(i,:));
        W_element = W(num_a1,:);
        WKW_free = WKW_free + W_element'*xpassive(i)'*KE(logcal_element,logcal_element)*W_element;
    end
    %% DCG 有限元分析
    K_gpu = gpuArray(K(freedofs,freedofs));
    W_gpu = gpuArray(W(freedofs,:));
    F1_gpu = gpuArray(F(freedofs,1));
    F2_gpu = gpuArray(F(freedofs,2));
    U1 = gpuArray(U1);
    U2 = gpuArray(U2);
    WKW_1 = WKW_free^-1;
    [U1(freedofs,:)]=DCG_WKW_1(K_gpu,W_gpu,F1_gpu,WKW_1);
    [U2(freedofs,:)]=DCG_WKW_1(K_gpu,W_gpu,F2_gpu,WKW_1);
    U1 = gather(U1);
    U2 = gather(U2);
    %% 更新设计变量
    [T1] = ComputeT(U1,x(:),edofMat,E0,nu); % 工况1拓扑导数
    [T2] = ComputeT(U2,x(:),edofMat,E0,nu); % 工况2拓扑导数
    T=0.5*T1+0.5*T2;                        % 拓扑导数加权
    T(x==0|passiveEleInd==1)=nan;
    xnew=x;
    if sum(sum(sum(x)))>=NumSolid*volfrac
        [~,I] = sort(T);
        IndexRemove=I(1:Tpercent);   % 确定去除单元
        xnew(IndexRemove)=0;         % 更新设计变量
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    ce = reshape(sum((U1(edofMat)*KE).*U1(edofMat),2),[nely,nelx,nelz])+reshape(sum((U2(edofMat)*KE).*U2(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum(x.*ce)));%结构柔顺度
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,sum(x(:))/NumSolid,change);% 输出当前迭代步结果
    %% 画图
    clf;
    xPhys = reshape((H*x(:))./Hs,nely,nelx,nelz);
    [f,v] = isosurface(xPhys,0.35);
    patch('Faces',f,'Vertices',v); alpha(0.5)
    axis equal; axis tight;  box on; view([45,50]);pause(1e-6);
    saveas(gcf,num2str(loop),'jpg');
    Record=[Record,x(:)];
end
save Record
clf; display_3D(x);
saveas(gcf,'simp','jpg');
end

%% DCG有限元分析子程序
function [x]=DCG_WKW_1(K,W,f,WKW_1)
x_1 = 0*ones(length(K),1);
r_1 = f-K*x_1;
x0 = x_1+W*(WKW_1*(W'*r_1));
r0 = f-K*x0;
mu0 = WKW_1*(W'*(K*r0));
p0 = r0 - W*mu0;
j=0;
while norm(r0)>1e-8
    Kp0 = K*p0;
    alpha0 = r0'*r0/(p0'*Kp0);
    x1 = x0 + alpha0*p0;
    r1 = r0 - alpha0*Kp0;
    bata0 = r1'*r1/(r0'*r0);
    mu1 = WKW_1*(W'*(K*r1));
    p1 = bata0*p0 + r1 -W*mu1;
    r0 = r1;
    x0 = x1;
    mu0 = mu1;
    p0 = p1;
    j=j+1;
end
x=x0;
end

%% 单元刚度矩阵子程序
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];
K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end

%% 画图子程序
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5) 
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight;  box on; view([30,30]); hold off;pause(1e-6);
end

%% 拓扑导数计算子程序
function [T] = ComputeT(U,ElmH,edofMat,E0,nu)
T = zeros(size(ElmH,1),1);
gradN =0.25 * [  1 -1 -1  1  1 -1 -1  1;
    1  1 -1 -1  1  1 -1 -1;
    -1 -1  1  1  1  1 -1 -1];
uGrad = gradN*U(edofMat(:,1:3:end))';
vGrad = gradN*U(edofMat(:,2:3:end))';
wGrad = gradN*U(edofMat(:,3:3:end))';
strains = [uGrad(1,:); vGrad(2,:); wGrad(3,:);uGrad(2,:)+vGrad(1,:);vGrad(3,:)+wGrad(2,:);uGrad(3,:)+wGrad(1,:)];
D0 = 1/(1+nu)/(1-2*nu)*[1-nu  nu   nu  0 0 0;nu  1-nu  nu  0 0 0;nu   nu  1-nu 0 0 0;
    0 0 0  (1-2*nu)/2 0 0;0 0 0  0 (1-2*nu)/2 0;0 0 0  0 0 (1-2*nu)/2;];
stresses = D0*E0*strains;
strains = reshape(strains,size(strains,1),1,size(strains,2));
stresses = reshape(stresses,size(stresses,1),1,size(stresses,2));
stressTensor = [stresses(1,:,:) stresses(4,:,:) stresses(6,:,:);
    stresses(4,:,:) stresses(2,:,:) stresses(5,:,:);
    stresses(6,:,:) stresses(5,:,:) stresses(3,:,:)];
strainTensor = [strains(1,:,:)   strains(4,:,:)/2 strains(6,:,:)/2;
    strains(4,:,:)/2 strains(2,:,:)   strains(5,:,:)/2;
    strains(6,:,:)/2 strains(5,:,:)/2 strains(3,:,:)];
for oneElm = 1:size(ElmH,1)
    alpha = ElmH(oneElm);
    if ElmH(oneElm) >= 0.5
        T(oneElm) =4/(1+nu)*sum(sum(stressTensor(:,:,oneElm).*strainTensor(:,:,oneElm)))- ...
            (1-3*nu)/(1-nu^2)*trace(stressTensor(:,:,oneElm))*trace(strainTensor(:,:,oneElm))*alpha;
    else
        T(oneElm) = 0;
    end
end
end
