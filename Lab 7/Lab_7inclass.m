clear; clc;
N = 20; nu=0.01;
x = linspace(0,1,N+1);
dx = x(2)-x(1);
a = (2*nu+dx)*ones(N+1,1);
b = -4*nu*ones(N+1,1);
c = (2*nu-dx)*ones(N+1,1);
d = zeros(N+1,1);
a(1) = 0; b(1) = 1; c(1) = 0; d(1) = 0;
a(end) = 0; b(end) = 1; c(end) = 0; d(end) = 1;

u = TDMA(a,b,c,d);
xx = linspace(0,1,1001);
plot(x,u,'-o',xx,(exp(xx/nu)-1)/(exp(1/nu)-1))