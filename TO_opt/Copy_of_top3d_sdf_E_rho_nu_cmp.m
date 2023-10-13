function [xPhys, yPhys, zPhys, output] = top3d_sdf_E_rho_nu_cmp(nelx,nely,nelz,volfrac,rmin)
clc;clear;
nelx=30;
nely=10;
nelz=1;
volfrac=0.45;
rmin=3;
% 用户定义循环参数
maxloop = inf;    % 最大迭代数
tolx = 0.00001;      % 终止条件
displayflag = 0;  % 结构显示
% 初始材料参数
E0 = 41.7;         % 初始化杨氏模量
nu0 = 0.2778;       % Poisson's ratio
rho0 = 0.5014;      %  density

%% add passive element 
passive = zeros(nely,nelx,nelz);
n = 0;
p = 0;
for k = 1:nelz
    for i = 1:nelx
        for j = 1:nely
            n = n+1;
            if sqrt((k-nelz)^2+(j-nely/2)^2+(i-nelx/3)^2)< nely/3
                passive(j,i,k) = 1;
                p = p+1;
                PassNum(p) = n;
            end
        end
    end
end
%有效单元数
ActEle = nelx*nely*nelz - size(PassNum,2);

% 载荷条件
[il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % 网格化
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % 节点ID
loaddof = 3*loadnid(:) - 1;                             % 自由度
% 固定自由度
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                 
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); 
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2];
% 准备有限元分析
nele = nelx*nely*nelz; % 单元数
ndof = 3*(nelx+1)*(nely+1)*(nelz+1); % 总自由度数
F = sparse(loaddof,1,-100,ndof,1); % 载荷条件
U = zeros(ndof,1); % 初始化位移向量
freedofs = setdiff(1:ndof,fixeddof); 
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1); 
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);

SimEdofMat=RemoveEdofMat(edofMat,PassNum);
SimIK = reshape(kron(SimEdofMat,ones(24,1))',24*24*ActEle,1);
SimJK = reshape(kron(SimEdofMat,ones(1,24))',24*24*ActEle,1);




% 准备灵敏度过滤
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
                        if passive(j2,i2,k2) == 1 
                            sH(k) = 0;
                        end
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH); 
Hs = sum(H,2); 

%% 删除H中的多余行和列
SimH=RemoveH(H,PassNum);
SimHs = sum(SimH,2);

% 去除消极单元的初始参数
X = [];
x = repmat(E0,[1,nelx*nely*nelz-size(PassNum,2)]);
y = repmat(nu0,[1,nelx*nely*nelz-size(PassNum,2)]);
z = repmat(rho0,[1,nelx*nely*nelz-size(PassNum,2)]);
X(1,:)=x;
X(2,:)=y;
X(3,:)=z;

x=reX(x,PassNum); xPhys = reshape(x,[nely,nelx,nelz]);
y=reX(y,PassNum); yPhys = reshape(y,[nely,nelx,nelz]);
z=reX(z,PassNum); zPhys = reshape(z,[nely,nelx,nelz]);
A = [];
B = [];
Aeq = [];
Beq = [];
LB = [0.01*ones(1,size(X,2)); 0.01*ones(1,size(X,2)); 0.01*ones(1,size(X,2))];
UB = [128*ones(1,size(X,2));  0.328*ones(1,size(X,2)); ones(1,size(X,2))];
options = optimoptions(@fmincon,'Display','iter','Algorithm','interior-point','StepTolerance',tolx,'MaxIterations',maxloop,'MaxFunctionEvaluations',inf,'ConstraintTolerance',1e-10,'HessianApproximation','bfgs',...
                       'OutputFcn',@(X,optimValues,state) myOutputFcn(X,optimValues,state,displayflag),'PlotFcn',@optimplotfval);


function f = myObjFcn(X)
    x = X(1,:);
    y = X(2,:);
    z = X(3,:);
    SimxPhys = (SimH*x(:))./SimHs;
    SimyPhys = (SimH*y(:))./SimHs;
    SimzPhys = (SimH*z(:))./SimHs;
    xPhysf = SimxPhys;
    yPhysf = SimyPhys;
    zPhysf = SimzPhys; 
    sK = [];
    for i = 1:size(yPhysf,1)
        KE = lk_H8(yPhysf(i));
        sK = horzcat(sK, KE(:)*xPhysf(i));
    end
    % 有限元分析
    sK = reshape(sK,24*24*ActEle,1);
    K = sparse(SimIK,SimJK,sK); K = (K+K')/2;A = full(K);
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % 目标函数评估
    UE = U(edofMat);
    ce = [];
    for j = 1:nely*nelx*nelz
        ce = [ce; (UE(j,:)*lk_H8(yPhysf(i))).*UE(j,:)];
    end
    ce = reshape(sum(ce,2),[nely,nelx,nelz]);
    c = sum(sum(sum(xPhys.*ce)));
    f = c + 0*(sum(zPhysf) - volfrac*nele)^2;
end 

function [cneq, ceq, gradc, gradceq] = myConstrFcn(X)
    x = X(:,1:nelx,:);
    y = X(:,nelx+1:2*nelx,:);
    z = X(:,2*nelx+1:3*nelx,:);
    xPhys(:) = (SimH*x(:))./SimHs;
    yPhys(:) = (SimH*y(:))./SimHs;
    zPhys(:) = (SimH*z(:))./SimHs;   
    % 非线性约束
    % 加载材料属性
    mat_prp = load('..\data\mat_prp.mat');
    mat_prp = mat_prp.mat_prp;
    % 加载点空间均值
    rr_ave = load('..\data\rr_ave_3d.mat');
    rr_ave = rr_ave.rr_ave;
    % 加载平均位置
    p_bar = load('..\data\p_bar_3d.mat');
    p_bar = p_bar.ans;
    % 归一化处理
    xs = mat_prp(:,1); xs_max = max(xs);
    % 
    p_ne = zeros(nelx*nely*nelz,3);
    for i = 1:nelx*nely*nelz
        kid = dsearchn(p_bar, [xPhys(i)/xs_max, yPhys(i), zPhys(i)]);
        p_ne(i,:) = p_bar(kid, :);
    end
    xp = p_ne(:,1);yp = p_ne(:,2);zp = p_ne(:,3);
    cneq = [];
    for i = 1:nelx*nely*nelz
        cneq1 = 1000*(((xPhys(i)/xs_max-xp(i))^2 + (yPhys(i)-yp(i))^2 + (zPhys(i)-zp(i))^2)^(0.5) - 8.0*rr_ave);
        cneq = [cneq cneq1];
    end
    cneq2 = sum(zPhys(:)) - volfrac*nele;
    cneq = [cneq cneq2];
    gradc = [];
    ceq     = [];
    gradceq = [];
end 

function stop = myOutputFcn(X,optimValues,state,displayflag)
    stop = false;
    switch state
        case 'iter'
            x = X(:,1:nelx,:);
            y = X(:,nelx+1:2*nelx,:);
            z = X(:,2*nelx+1:3*nelx,:);
            xPhys = reshape(x, nely, nelx, nelz);
            yPhys = reshape(y, nely, nelx, nelz);
            zPhys = reshape(z, nely, nelx, nelz);
            %% 输出结果
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',optimValues.iteration,optimValues.fval, ...
                mean(zPhys(:)),optimValues.stepsize);
            %% 绘制结果图
            if displayflag, figure(10); clf; display_3D(zPhys); end
            title([' It.:',sprintf('%5i',optimValues.iteration),...
                ' Obj. = ',sprintf('%11.4f',optimValues.fval),...
                ' ch.:',sprintf('%7.3f',optimValues.stepsize)]);

        case 'init'
            if displayflag
                figure(10)
            end
        case 'done'
            figure(10); clf; display_3D(zPhys);
            save compdata
        otherwise
    end 
end 

output = fmincon(@(X)myObjFcn(X), X, A, B, Aeq, Beq, LB, UB, @(X)myConstrFcn(X), options);
end

% === 生成单元刚度矩阵 ===
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
% === 3D结果显示程序 ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;           
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'color','w','Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.10)  
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'LineStyle','none','FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))],'FaceAlpha',.2);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end


%% === 插入消极单元，还原源矩阵 ===
function CompX=reX(SimX,PassNum)
x=SimX;
    for i = 1:size(PassNum,2)
        temp=x(PassNum(i):end);
        x(PassNum(i)) = nan;
        x(PassNum(i)+1:end+1) = temp;
    end
CompX=x;
end


%% === 删除H中的多余行和列 ===
function SimH=RemoveH(H,PassNum)
CompH = H;
    for i = 1:size(PassNum,2)
        temp1=CompH(PassNum(i)+1:end,:);
        temp2=CompH(:,PassNum(i)+1:end);
        CompH(PassNum(i):end-1,:) = temp1;
        CompH(:,PassNum(i):end-1) = temp2;
        CompH(end,:) = [];
        CompH(:,end) = [];
    end
    SimH = CompH;
end


%% === 删除edofMat中的多余行 ===
function SimEdofMat=RemoveEdofMat(edofMat,PassNum)
CompEdofMat = edofMat;
SimEdofMat = CompEdofMat;
    for i = 1:size(PassNum,2)
        temp=CompEdofMat(PassNum(i)+1:end,:);
        SimEdofMat(PassNum(i)-i+1:end-1-i+1,:) = temp;
    end
SimEdofMat(end-size(PassNum,2)+1:end,:) = [];
end