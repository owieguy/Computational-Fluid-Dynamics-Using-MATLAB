clear; clc;
N=101; x=linspace(0,1,N+1); dx=x(2)-x(1);
m=zeros(N+1,1); s=zeros(N+1,1);
for i = 1:N
    if x(i) >= 0.3 && x(i) <= 0.7
        m(i)=1;s(i)=1;
    end 
end 

u=ones(N+1,1);p=zeros(N+1,1);
a=ones(N+1,1); a(1)=0; a(N+1)=0;
b=-2*ones(N+1,1); b(1)=-1; b(N+1)=1;
c=ones(N+1,1); c(1)=1; c(N+1)=0;
d=zeros(N+1,1); d(1)=0; d(N+1)=0;
for n=1:30000 % Time marching
    dt=min(0.5*dx^2,min(2./u.^2));
    for i = 2:N
        d(i)=-m(i)^2*dx^2.....
            -u(i)*(u(i-1)-2*u(i)+u(i+1))....
            +m(i-1)-2*m(i)+m(i+1).....
            +dx/2*(s(i+1)-s(i-1));
    end
    p=Step_05_Poissons_Equation(a,b,c,d);
    un=u;
    u(1)=1;
    for i=2:N
        u(i)=un(i)....
            -(un(i)*dt/2/dx)*(un(i+1)-un(i-1))...
            +(dt/dx/dx)*(un(i-1)-2*un(i)+un(i+1))....
            -(dt/2/dx)*(p(i+1)-p(i-1))...
            +dt*s(i);
    end
    u(N+1)=u(N);
    if ceil(n/100)==n/100
        plot(x,u,'o'),pause(0.02);
        plot(x,p,'o'),
    end 
end

