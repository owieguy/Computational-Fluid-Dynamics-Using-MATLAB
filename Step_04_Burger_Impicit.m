clear; clc;
nu=1; N=100; 
x=linspace(-2,4,N+1); dx=0.06;
u=ones(N+1,1);
nt=101;
for i=1:N+1
    if x(i)>0.5 && x(i)<1
        u(i)=2;
    end
end
dt=min(0.5*dx^2/nu,min(2*nu./u.^2));

for n=1:nt
    un=u;
    err=1;
    while err>1e-5
        um=u;
        a = -u*dt/2/dx-nu*dt/dx^2*ones(1,N+1);
        b = 1+2*nu*dt/dx^2*ones(1,N+1);
        c = u*dt/2/dx-nu*dt/dx^2*ones(1,N+1);
        d = un;
        a(1)=0; b(1)=1; c(1)=0; d(1)=1;
        a(N+1)=0; b(N+1)=1; c(N+1)=0; d(N+1)=1;
        u = TDMA(a,b,c,d); %Um+1

        err = norm(u-um)
    end 
     plot(x,u),
     axis([-2 4 1 2]),
     pause(0.2);
end