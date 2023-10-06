clear; clc;
N = 20; c=1; new=0.01;
x = linspace(0,2,N+1);
tf = 3; dx = x(2)-x(1);dt = 0.25*min((dx/c),(dx^2/new));
M=tf/dt; sigma=c*dt/dx; small=1e-8;
%+++++++++++++++++++++%
u = ones(N+1,1);
%+++++++++++++++++++%

for n =1:M
    un = u;
    for i= 2:N
        r(i)=(u(i)-u(i-1))/(u(i+1)-u(i)+small);
    end
    r(1)= r(2);
    r(N+1)= r(N);
    phi=max(max(0,min(2*r,1)),min(r,2)); %Superbee limiter
    for i=2:N
        u(i) = un(i) - sigma*(un(i)-un(i-1))...
            -sigma/2*(1-sigma)*(phi(i)/(r(i)+small)-phi(i-1))*(u(i)-un(i-1))+(new*dt/dx^2)*(u(i+1)-2*u(i)+u(i-1));
    end
    u(1) = 0;
    u(N+1) = 1;
    plot(x,u,'-o')                                                                                                                                                                                                                                                                                                                                                                                                                                             
    pause(0.2)
end