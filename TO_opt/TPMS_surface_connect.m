function []=TPMS_surface_connect(alpha1_now,alpha2_now,alpha3_now,alpha1_old,alpha2_old,alpha3_old,t1,t2,t3,Xc,Yc,Zc)
    % alpha1 = 0.411;
    % alpha2 = 0.144;
    % alpha3 = 0.445;
    % t1 = -0.120;
    % t2 = -0.349;
    % t3 = -0.359;
    
    x = [-0.5:0.01:0.5]*2*pi;
    y = [-0.5:0.01:0.5]*2*pi;
    z = [-0.5:0.01:0.5]*2*pi;
    [X,Y,Z] = meshgrid(x,y,z);
    
    F1_now = cos(X) + cos(Y) + cos(Z)+ t1;
    F2_now = cos(X).*cos(Y).*cos(Z) - sin(X).*sin(Y).*sin(Z) + t2;
    F3_now = 8*cos(X).*cos(Y).*cos(Z)+cos(2*X).*cos(2*Y).*cos(2*Z)-(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))+t3;

    n1 = alpha1_now+alpha1_old * X;
    n2 = alpha2_now+alpha2_old * X;
    n3 = alpha3_now+alpha3_old * X;

    Fmerge = n1.*4.*F1_now+n2.*4.*F2_now+n3.*F3_now;

    
    % fun = @(X,Y,Z)0.576*(cos(X)+cos(Y)+cos(Z)+0.146)+0.364*(cos(X).*cos(Y).*cos(Z)-sin(X).*sin(Y).*sin(Z)+-0.0547)+0.0597*(8*cos(X).*cos(Y).*cos(Z)+cos(2*X).*cos(2*Y).*cos(2*Z)-(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))+-0.213);
    
    figure(1001)
    isosurface(X+2*pi*Xc,Y-2*pi*Yc,Z-2*pi*Zc,Fmerge,1e-4)
    axis equal;
end