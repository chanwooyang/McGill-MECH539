function [ u ] = Lax( dx, dt )
%INPUT
%dt: Discretization size of time domain
%dx: Discretization size of spatial domain

%OUTPUT
%u: Array of values of wave equation at each grid point


x=0;    %Spatial domain boundary (start)
X=40;   %Spatial domain boundary (end)
t=0;    %Time domain boundary (start)
T=10;   %Time domain boundary (end)
C=1/2;  %C = wave speed
v=C*dt/dx;    

% Initialize the exact Dirichlet boundary conditions
for n=1:((T/dt)+1)
    u(1,n)=0;
    u(((X/dx)+1),n)=1;
end

% Initialize the initial condition
for j=1:((X/dx)+1)
    u(j,1)=(1/2)*(1+tanh(250*(x-20)));
    x=x+dx;
end

%Lax Method
for n=1:(T/dt)
    for j=2:(X/dx)
        u(j,n+1)=(1-v)*(1/2)*u(j+1,n)+(1+v)*(1/2)*(u(j-1,n));
    end
end

end