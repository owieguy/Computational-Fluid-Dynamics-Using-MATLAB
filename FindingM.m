clear; clc; 
R = 287; cv = 718; cp = 1005; k = 1.4; %Parameters of Air
N = 60; x = linspace(0,3,N+1); dx = x(2) - x(1);
A = 1 + 2.2*(x-1.5).^2; 
p0 = 1e5; pe = 0.6784e5;
p = linspace(p0,pe,N+1);
rho = ones(1,N+1);
V = zeros(1,N+1);
T = p./(rho*R);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
U = findU(rho,A,T,V); F = findF(rho,A,T,V); S = findS(p,x);
Ubar=U; Ubarbar=U;
err = 1; i = 2 : N;
while err > 1e-3
    Vn=V;
    dt = dx/max(abs(V+sqrt(k*R*T)));
    Un = U;
    %first step
    Ubar(:,i) = Un(:,i).....
        -(dt/dx)*(F(:,i)-F(:,i-1)).......
        +dt*S(:,i);
    Ubar(1,1) = U(1,1); %Specify U(1) at entrance
    Ubar(2,1) = 2*Ubar(2,2) - Ubar(2,3); %Extrapolate U(2)
    Ubar(3,1) = U(3,1); %Specify U(3)

    Ubar(1,N+1) = 2*Ubar(1,N) - Ubar(1,N-1); %Extrapolate U(1)
    Ubar(2,N+1) = 2*Ubar(2,N) - Ubar(2,N-1); %Extrapolate U(2)
    Ubar(3,N+1) = U(3,N+1); %Specify U(3) at the exit

    [rho,T,V,p] = findVar(Ubar,A);
    Fbar = findF(rho,A,T,V); Sbar = findS(p,x);

    %Second Step
    Ubarbar(:,i) = Un(:,i).....
        -(dt/dx)*(Fbar(:,i+1)-Fbar(:,i)).......
        +dt*Sbar(:,i);
    Ubarbar(1,1) = U(1,1); %Specify U(1) at entrance
    Ubarbar(2,1) = 2*Ubarbar(2,2) - Ubarbar(2,3); %Extrapolate U(2)
    Ubarbar(3,1) = U(3,1); %Specify U(3)

    Ubarbar(1,N+1) = 2*Ubarbar(1,N) - Ubarbar(1,N-1); %Extrapolate U(1)
    Ubarbar(2,N+1) = 2*Ubarbar(2,N) - Ubarbar(2,N-1); %Extrapolate U(2)
    Ubarbar(3,N+1) = U(3,N+1); %Specify U(3) at the exit
    %Third step
    U = 0.5*(Ubar+Ubarbar);
    [rho,T,V,] = findVar(U,A);
    F = findF(rho,A,T,V); S = findS(p,x);
    err=norm(V-Vn);
    plot(x,V./sqrt(k*R*T),'-o'); pause(0.1);
    %plot(x,V./sqrt(k*R*T),'o',x,FindExact(x,A)); pause(0.1);
end
plot(x,V./sqrt(k*R*T),'o');
function U = findU(rho,A,T,V)
R = 287; cv = 718; cp = 1005; k = 1.4;
U(1,:) = rho.*A;
U(2,:) = rho.*A.*V;
U(3,:) = rho.*A.*(cv*T+0.5*V.^2);
end
function F = findF(rho,A,T,V)
R = 287; cv = 718; cp = 1005; k = 1.4;
p=rho.*R.*T;
F(1,:) = rho.*A.*V;
F(2,:) = rho.*A.*V.^2+p.*A;
F(3,:) = rho.*A.*V.*(cp*T+0.5*V.^2);
end
function S = findS(p,x)
R = 287; cv = 718; cp = 1005; k = 1.4;
S(1,:) = zeros(size(p));
S(2,:) = p.*4.4.*(x-1.5);
S(3,:) = zeros(size(p));
end
function [rho,T,V,p] = findVar(U,A)
R = 287; cv = 718; cp = 1005; k = 1.4;
rho = U(1,:)./A;
V = U(2,:)./(U(1,:));
T = (U(3,:)./U(1,:)-0.5*V.^2)/cv;
p = rho*R.*T;
end