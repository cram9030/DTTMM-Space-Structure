function [A,B,D,E] = Houbolt(dt,x,dx,ddx,x2,ddx2,x3,i)

A = 2/(dt^2);
if i>=3
    D = 11/(6*dt);
    B = (-1/dt^2)*(5*x-4*x2+x3);
    E = (-1/(6*dt))*(18*x-9*x2+2*x3);
elseif i==2
    D = 11/(6*dt);
    B = (-1/dt^2)*(4*x-2*x2+2*dt^2*ddx2);
    E = (-1/(6*dt))*(16*x-5*x2+dt^2*ddx2);
elseif i==1
    D = 3/dt;
    B = (-2/dt^2)*(3*x+3*dt*dx+dt^2*ddx);
    E = (-1/(2*dt))*(6*x+4*dt*dx+dt^2*ddx);
end