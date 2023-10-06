N=101;
x = linspace(0,1,N+1);
dx=x(2)-x(1);
a=ones(1,N+1); a(1)=0; a(N+1)=-1;
b=-(2+dx^2)*ones(1,N+1); b(1)=1; b(N+1)=1;
c=ones(1,N+1); c(1)=0; c(N+1)=0;
d=zeros(1,N+1); d(N+1)=dx;
u=TDMA(a,b,c,d);
plot  (x,u,'o')
T=u(N+1)
function u=TDMA(a,b,c,d)
% a, b, c, d are vectors, which contain
% [a1, a2, a3,....., aN]
N=length(a); u=zeros(N,1); 
q=u; p=u; 
q(1) = d(1)/b(1); p(1) = c(1)/b(1);
for i = 2:N
    q(i) = (d(i)-a(i)*q(i-1))/(b(i)-a(i)*p(i-1));
    p(i) = c(i)/(b(i)-a(i)*p(i-1));
end
u(N) = q(N);
for i = N-1:-1:1
    u(i)=q(i)-p(i)*u(i+1);
end 
end