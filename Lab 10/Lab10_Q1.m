clear; clc;
nu = 0.02; delta = 0.08;
N1 = 24; x1 = linspace(0,0.92,N1+1);
N2 = 6; x2 = linspace(0.92,1,N2+1);
N = N1+N2;
x = [x1,x2(2:end)]';

Delta_x = diff(x);
Delta_x = [Delta_x(1);Delta_x; Delta_x(end)];

xc = 0.5*(x(1:end-1)+x(2:end));
xc = [xc(1)-Delta_x(1);xc;xc(end)+Delta_x(end)];
dx = diff(xc);
a = zeros(N+2,1); b=a; c=a;d=a;Face1=a; Face2=a;
u=ones(N+2,1);


err=1;
while err>0.000001
    un = u;
    i = 2:N+1;
    
    Face1(i)=(u(i)+u(i+1))/4;
    Face2(i)=(u(i-1)+u(i))/4;
    a(i) = nu./dx(i-1)+Face2(i);
    b(i) = -(nu*(1./dx(i-1)+1./dx(i)))....
    -Face1(i)+Face2(i);
    c(i) = nu./dx(i)-Face1(i);

    a(1) = 0; b(1) = 0.5; c(1) = 0.5; d(1) = 1;
    a(N+2) = 0.5; b(N+2) = 0.5; c(N+2) = 0; d(N+2) = 0;
    

    u = TDMA(a,b,c,d);
    err = max(abs(un-u));
end
%(x(i-1)-x(i))./(x(i)-x(i+1))
plot(xc,u,'-o')