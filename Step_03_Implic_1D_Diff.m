clear; clc; 
N=200; nu=1; 
x=linspace(0,1,N+1); dx=x(2)-x(1);
dt=0.5*dx^2/nu;  tf=0.05;
M=tf/dt;

u=zeros(N+1,1); u(N+1)=1;

r=nu*dt/dx^2

a=-r*ones(1,N+1); b=(1+2*r)*ones(1,N+1);
c=-r*ones(1,N+1); d=zeros(1,N+1);
a(1)=0; b(1)=1; c(1)=0; d(1)=1;
a(N+1)=0; b(N+1)=1; c(N+1)=0; d(N+1)=1;

for n=1:M
    un=u;
    for i = 2:N
        u(i) = un(i)+(nu*dx/dx^2)*(un(i+1)-2*un(i)+un(i-1));
    end
    for i = 1:N+1
        ue(i)=u_exact(x(i),n*dt);
    end

    
    plot(x,u,'o',x,ue)
    axis([-0 1 0 1])
    pause(0.2)
end

function res = u_exact(x,t)
m=1:100;
res=x-2*sum(exp(-(m*pi).^2*t).*sin(m*pi*(1-x))./(m*pi));
end