clear; clc;
nu=1; N=101; x=linspace(-2,4,N+1); dx=0.06;
u=ones(N+1,1);

for i=1:N+1
    if x(i)>0.5 && x(i)<1
        u(i)=2;
    end
end
%dt=min(0.5*dx^2/nu,min(2*nu./u.^2))
dt=0.002
for n=1:500
    un=u;
    for i=2:N
        u(i)=un(i)-un(i)*dt/2/dx*(un(i+1)-un(i-1))...
            +nu*dt/dx^2*(un(i+1)-2*un(i)+un(i-1));
    end
     u(1)=un(1)-un(1)*dt/2/dx*(un(2)-un(N))...
            +nu*dt/dx^2*(un(2)-2*un(1)+un(N));
     u(N+1)=un(N+1)-un(N+1)*dt/2/dx*(un(2)-un(N))...
            +nu*dt/dx^2*(un(2)-2*un(N+1)+un(N));
     plot(x,u),
     axis([-2 4 1 2]),pause(0.2);
end