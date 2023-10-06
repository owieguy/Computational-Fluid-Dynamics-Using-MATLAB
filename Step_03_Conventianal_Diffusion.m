clear; clc;
nu = 0.001; delta = 6*nu;
N1 = 4; x1 = linspace(0,1-delta,N1+1);
N2 = 4; x2 = linspace(1-delta,1,N2+1);
N = N1+N2;
x = [x1,x2(2:end)];
Delta_x = diff(x);
Delta_x = [Delta_x(1),Delta_x,Delta_x(end)];

xc = 0.5*(x(1:end-1)+x(2:end));
xc = [xc(1)-Delta_x(1),xc,xc(end)+Delta_x(end)];
dx = diff(xc);
a = zeros(N+2,1); b=a; c=a;d=a;
i = 2:N+1;
a(i) = 2*nu./dx(i-1)+1;
b(i) = -2*nu*(1./dx(i-1)+1./dx(i));
c(i) = 2*nu./dx(i)-1;

a(1) = 0; b(1) = 0.5; c(1) = 0.5; d(1) = 0;
a(N+2) = 0.5; b(N+2) = 0.5; c(N+2) = 0; d(N+2) = 1;

u = TDMA(a,b,c,d);
xx = linspace(0,1,1001);
plot(xc,u,'-o',xx,(exp(xx/nu)-1)/(exp(1/nu)-1))
error = max(abs(u(2:end-1)'-(exp(xc(2:end-1)/nu)-1)/(exp(1/nu)-1)))