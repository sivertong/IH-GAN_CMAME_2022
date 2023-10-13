function []=remakeK()
    clc;clear;
% nelx=30;
% nely=10;
% nelz=1;
% volfrac=0.5;
% rmin=3;

nelx=10;
nely=10;
nelz=1;
volfrac=0.5;
rmin=3;
% USER-DEFINED LOOP PARAMETERS
maxloop = inf;    % Maximum number of iterations
tolx = 0.00001;      % Terminarion criterion
displayflag = 0;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E0 = 41.7;         % Initial Young's modulus
nu0 = 0.2778;       % Initial Poisson's ratio
rho0 = 0.5014;      % Initial density
% E0 = 30.17;         % Initial Young's modulus
% nu0 = 0.3075;       % Initial Poisson's ratio
% rho0 = 0.4513;      % Initial density
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz; % Number of elements
ndof = 3*(nelx+1)*(nely+1)*(nelz+1); % Total number of degree of freedom
F = sparse(loaddof,1,-100,ndof,1); % Apply loading condition
U = zeros(ndof,1); % Initialize displacement
freedofs = setdiff(1:ndof,fixeddof); % Indices indicating unconstrained DOFs (excluding the fixed DOFs)
% KE = lk_H8(nu); % Element stiffness matrix (24 by 24)
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1); % Node indices for each element following the local node order (nele by 24)

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
                elemNum(p) = n;
            end
        end
    end
end



iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1); % Indices for global stiffness matrix




end

