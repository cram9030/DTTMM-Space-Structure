%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,x_0,dx_0,ddx,x2,ddx2,x3,iter,fn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DTTMMHou2015_1_30: Single Discrete Time Transfer Method solving time step	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves a sinlge time step for a Galerkin Finite element
% discrete time transfer matrix method using the Hubolt integration scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:	DTFETMMHou_Loop_SS_2015_1_30									
% Calls:	MATLAB 5.2 std fcns
% Inputs:   dt - Time step
%           E -   Modulus of elasticity (Pa)
%           etaE = Damping coefficient (Pa/s)
%           I - Area moment of inertia (m^4)
%           rho - Density (kg/m^3)
%           CA - Cross sectional area (m^2)
%           L - Length array (m)
%           x_0 - Previous time step position array
%           dx_0 - Previous time step velocity array
%           ddx - Previous time step acceleration array
%           x2 - Two time steps ago position array
%           dx2 - Two time steps ago velocity array
%           ddx2 - Two time steps ago acceleration array
%           x3 - Three time steps ago position array
%           iter - Current iteration to know what version of Houbolt
%           integration to use
%           fn - Array of forces
% Outputs:  x - current state vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initalize Matrices
SegNum = length(x_0);
TL = eye(5);
TR = eye(5);
Q = eye(5);
G = eye(5);
x = zeros(4,length(x_0)+1);
x_0 = [zeros(2,1) x_0 zeros(2,1)];
dx_0 = [zeros(2,1) dx_0 zeros(2,1)];
ddx = [zeros(2,1) ddx zeros(2,1)];
x2 = [zeros(2,1) x2 zeros(2,1)];
ddx2 = [zeros(2,1) ddx2 zeros(2,1)];
x3 = [zeros(2,1) x3 zeros(2,1)];
fn = [zeros(2,1) fn zeros(2,1)];
P = zeros(5,5,length(x_0)-1);
F = zeros(5,5,length(x_0)-1);

%Populate P and F matrices
for i = 2:SegNum+1
    %Create stiffness sub matrices
    K11 = E*I/L(i-1)^3*[12 6;6 4];
    K12 = E*I/L(i-1)^3*[-12 6;-6 2];
    K21 = E*I/L(i-1)^3*[-12 -6;6 2];
    K22 = E*I/L(i-1)^3*[12 -6;-6 4];

    %Create damping sub matrices
    C11 = etaE*I/L(i-1)^3*[12 6;6 4];
    C12 = etaE*I/L(i-1)^3*[-12 6;-6 2];
    C21 = etaE*I/L(i-1)^3*[-12 -6;6 2];
    C22 = etaE*I/L(i-1)^3*[12 -6;-6 4];

    %Create mass sub matrices
    M11 = rho*CA/(420*L(i-1))*[156 22;22 4];
    M12 = rho*CA/(420*L(i-1))*[54 -13;13 -3];
    M21 = rho*CA/(420*L(i-1))*[54 13;-13 -3];
    M22 = rho*CA/(420*L(i-1))*[156 -22;-22 4];
    Mnn = M11+M22;

    %Intial integration
    [An,Bn,Dn,En] = Houbolt(dt,x_0,dx_0,ddx,x2,ddx2,x3,iter);

    %Create constants that will be used for population of P and F matrices
    Finv = pinv(-C12*Dn-K12-M12*An);
    Fconst = C22*Dn+K22;

    BnL = Bn(:,i-1);
    EnL = En(:,i-1);
    Bnc = Bn(:,i);
    Enc = En(:,i);
    P(:,:,i-1) = [eye(2) zeros(2,3);Mnn*An eye(2) Mnn*Bnc-fn(:,i);zeros(1,4),1];
    F(:,:,i-1) = [Finv*(C11*Dn+K11),Finv,Finv*(C12*EnL+C12*Enc+M12*Bnc);...
        Fconst*Finv*(C11*Dn+K11)+C21*Dn+K21+M21*An,Fconst*Finv,C22*Enc+C21*EnL+M21*BnL+...
        Fconst*Finv*(C11*EnL+C12*Enc+M12*Bnc);...
        zeros(1,4),1];
end

%Create Q
Q = eye(5);
for i = 1:SegNum
        Q = P(:,:,i)*F(:,:,i)*Q;
end

%Solve for boundary conditions
x(3:4,1) = Q(3:4,3:4)\(-Q(3:4,5));
x(1:2,end) = Q(1:2,3:4)*x(3:4,1)+Q(1:2,5);

%Solve for the rest of the states
for i = 1:SegNum-1
    temp = P(:,:,i)*F(:,:,i);
    x(:,i+1) = temp(1:4,:)*[x(:,i);1];
end