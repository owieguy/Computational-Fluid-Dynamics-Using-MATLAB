clear; clc; clf;
L=1; k=50; pi=20; pe=0;
N = 400;
x = linspace(0,1,N+1); dx = x(2)-x(1);
A = linspace(0.5,0.1,N+1);
u = ones(N+1,1); %p = zeros(N+1,1);
a=ones(N+1,1).*(0.4/(2*dx)+A./dx^2); a(1)=0; a(N+1)=0;
b=ones(N+1,1).*(-2.*A./dx^2); b(1)=1; b(N+1)=1;
c=ones(N+1,1).*(-0.4/(2*dx)+A./dx^2); c(1)=0; c(N+1)=0;
d=zeros(N+1,1);

err = 1;
while err > 1e-12
%for j=1 : 100000
    dt = dx/max(k*u);
    un = u;
    for i = 1 : N+1
            u(i) = un(i)*(1-k*dt*un(i)); %This is u*
    end
    for i = 2 : N
        d(i) = (1/dt)*(A(i+1)*un(i+1)-A(i-1)*un(i-1))/(2*dx);
    end
    d(1)=20; d(N+1)=0;
    % now solve PPE
    p=TDMA(a,b,c,d);
    for i = 2 : N   
        u(i) = u(i) - dt*(p(i+1)-p(i-1))/2/dx;
    end
    u(1) = u(1) - dt*(p(2)-p(1))/dx;
    u(N+1) = u(N+1) - dt*(p(N+1)-p(N))/dx;
    err = abs((u(N+1)-un(N+1))/u(N+1));
end
plot(x,p);
u
u(N+1)