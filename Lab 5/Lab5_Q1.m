%clear; clc;
NM=200;
N=NM; M=NM; nu=2*10^-5;
x=linspace(0,0.1,N+1); dx=x(2)-x(1);
y=linspace(0,0.1,M+1); dy=y(2)-y(1);

[Y,X]=meshgrid(y,x);
u=ones(N+1,M+1);

for j=1:M+1
    for i=1:N+1
        if x(i)>=0.5 && x(i)<=1 && y(j)>=0.5 && y(j)<=1
            u(i,j)=2;
        end
    end 
end
dt08=0.5*min(dx/0.1,dy/0.1);
dt09=0.5/(nu*(1/dx^2+1/dy^2));
dt=min(dt08,dt09)*0.5;
error=10;
while abs(error) > 1E-10
    un=u;
    for j=2:M
        for i=2:N
            u(i,j)=un(i,j)....
                -(x(i)*dt/dx)*(un(i,j)-un(i-1,j))...
                +(y(j)*dt/dy)*(un(i,j+1)-un(i,j))...
                +(nu*dt/dx^2)*(un(i-1,j)-2*un(i,j)+un(i+1,j))...
                +(nu*dt/dy^2)*(un(i,j-1)-2*un(i,j)+un(i,j+1));
        end 
    end
    u(1,:)=u(2,:)+4000*dx; 
    u(:,1)=u(:,2)+4000*dy; 
    u(N+1,:)=u(N); 
    u(:,N+1)=10;
    error=u(1,1)-un(1,1);
end
Final=u(1,1)
%surf(X,Y,u)