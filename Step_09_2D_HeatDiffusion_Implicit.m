clear; clc;

N=80; M=80; nu=1;
x=linspace(-2,4,N+1); dx=x(2)-x(1);
y=linspace(-2,4,M+1); dy=y(2)-y(1);

[Y,X]=meshgrid(y,x);
u=ones(N+1,M+1);

for j=1:M+1
    for i=1:N+1
        if x(i)>=0.5 && x(i)<=1 && y(j)>=0.5 && y(j)<=1
            u(i,j)=2;
        end
    end 
end
dt=0.5/(nu*(1/dx^2+1/dy^2));
nt = ceil(0.5/dt);
rx = nu*dt/dx^2; ry = nu*dt/dy^2;
v = u;
for n=1:nt
    
    un=u; 

%Step 1 +++++++++++++++++++++++++++++++
    a = -rx*ones(N+1,1);
    b = 1+2*rx*ones(N+1,1);
    c = -rx*ones(N+1,1);
    a(1)=0; b(1)=1; c(1)=0; d(1)=1;
    a(N+1)=0; b(N+1)=1; c(N+1)=0; d(N+1)=1;
    for j = 2:M
        d(2:N) = un(2:N,j);
        v(:,j) = TDMA(a,b,c,d);
    end

%Step 2 +++++++++++++++++++++++++++++++
    a = -ry*ones(N+1,1);
    b = 1+2*ry*ones(N+1,1);
    c = -ry*ones(N+1,1);
    a(1)=0; b(1)=1; c(1)=0; d(1)=1;
    a(M+1)=0; b(M+1)=1; c(M+1)=0; d(M+1)=1;
    for i = 2:M
        d(2:M) = v(i,2:M);
        u(i,:) = TDMA(a,b,c,d);
    end

    surf(X,Y,u),
    axis([-2 4 -2 4 1 2]);
    pause(0.2);
end