function []=TPMS_surface(alpha1,alpha2,alpha3,t1,t2,t3,Xc,Yc,Zc)
% alpha1 = 0.411;
% alpha2 = 0.144;
% alpha3 = 0.445;
% t1 = -0.120;
% t2 = -0.349;
% t3 = -0.359;

x = [-0.5:0.1:0.5]*2*pi;
y = [-0.5:0.1:0.5]*2*pi;
z = [-0.5:0.1:0.5]*2*pi;
[X,Y,Z] = meshgrid(x,y,z);

F1 = cos(X) + cos(Y) + cos(Z)+ t1;
F2 = cos(X).*cos(Y).*cos(Z) - sin(X).*sin(Y).*sin(Z) + t2;
F3 = 8*cos(X).*cos(Y).*cos(Z)+cos(2*X).*cos(2*Y).*cos(2*Z)-(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))+t3;
Fmerge = alpha1*4*F1+alpha2*4*F2+alpha3*F3;
% fun = @(X,Y,Z)0.576*(cos(X)+cos(Y)+cos(Z)+0.146)+0.364*(cos(X).*cos(Y).*cos(Z)-sin(X).*sin(Y).*sin(Z)+-0.0547)+0.0597*(8*cos(X).*cos(Y).*cos(Z)+cos(2*X).*cos(2*Y).*cos(2*Z)-(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))+-0.213);

figure(1001)
isosurface(X+2*pi*Xc,Y-2*pi*Yc,Z-2*pi*Zc,Fmerge,1e-4)
axis equal;
end