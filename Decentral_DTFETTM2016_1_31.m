%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = Decentral_DTFETTM2016_1_31(dt,E,etaE,I,rho,CA,L,SegNum,node)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decentral_DTFETTM2016_1_31: Create decentralized model of space structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is called to create an decentralize model of a three node
% DTFETMM with input matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:										
% Calls:	MATLAB 5.2 std fcns
% Inputs:   dt - Time step
%           E - Modulus of elasticity
%           etaE - Material damping coefficent
%           I - Area moment of inertia of the model
%           rho - Density of the beam
%           CA - Cross sectional area
%           L - Array of section lengths
%           SegNum - Numbers of segments
%           node - Node to use as center node
% Outputs:  A - State matrix
%           B - Input matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Do not pass anything nodes within one node of an edge
%       SI units do not seem to work, pass values with kg-s-mm
%       Assume that the boundary conditions are clamped free

%Intialize matrices
A = zeros(18,30);
B = zeros(18,6);
int = 5;
TL = eye(5);
TR = eye(5);
Q = eye(5);
G = eye(5);
x = zeros(4,SegNum+1);
x_0 = zeros(2,SegNum+2);
dx_0 = zeros(2,SegNum+2);
ddx = zeros(2,SegNum+2);
x2 = zeros(2,SegNum+2);
ddx2 = zeros(2,SegNum+2);
x3 = zeros(2,SegNum+2);
P = zeros(5,5,length(x_0)-1);
F = zeros(5,5,length(x_0)-1);
H = zeros(5,5,length(x_0)-1);
J = zeros(5,5,length(x_0)-1);

%Generate P, F, H, and J
for i = 2:SegNum+1
    K11 = E*I/L(i-1)^3*[12 6;6 4];
    K12 = E*I/L(i-1)^3*[-12 6;-6 2];
    K21 = E*I/L(i-1)^3*[-12 -6;6 2];
    K22 = E*I/L(i-1)^3*[12 -6;-6 4];
    Knn = K11+K22;

    C11 = etaE*I/L(i-1)^3*[12 6;6 4];
    C12 = etaE*I/L(i-1)^3*[-12 6;-6 2];
    C21 = etaE*I/L(i-1)^3*[-12 -6;6 2];
    C22 = etaE*I/L(i-1)^3*[12 -6;-6 4];
    Cnn = C11+C22;

    M11 = rho*CA/(420*L(i-1))*[156 22;22 4];
    M12 = rho*CA/(420*L(i-1))*[54 -13;13 -3];
    M21 = rho*CA/(420*L(i-1))*[54 13;-13 -3];
    M22 = rho*CA/(420*L(i-1))*[156 -22;-22 4];
    Mnn = M11+M22;

    [An,Bn,Dn,En] = Houbolt(dt,x_0,dx_0,ddx,x2,ddx2,x3,int);
    Finv = (-C12*Dn-K12-M12*An)\eye(2);
    Hinv = (-C21*Dn-K21-M21*An)\eye(2);

    BnL = Bn(:,i-1);
    EnL = En(:,i-1);
    Bnc = Bn(:,i);
    Enc = En(:,i);
    Bnn = Bn(:,i+1);
    Enn = En(:,i+1);
    if i == 2
        Fconst = C11*Dn+K11;
        Hconst = C11*Dn+K11;
        P(:,:,i-1) = [eye(2) zeros(2,3);M11*An eye(2) M11*Bnc;zeros(1,4),1];
        F(:,:,i-1) = [Finv*(C11*Dn+K11),Finv,Finv*(C12*EnL+C12*Enc+M12*Bnc);...
            Fconst*Finv*(C11*Dn+K11)+C21*Dn+K21+M21*An,Fconst*Finv,Cnn*Enc+C21*EnL+M21*BnL+...
            Fconst*Finv*(C11*EnL+C12*Enc+M12*Bnc);...
            zeros(1,4),1];
        J(:,:,i-1) = [eye(2) zeros(2,3);-M11*An eye(2) -M11*Bnc;zeros(1,4),1];
        H(:,:,i-1) = [Hinv*(Cnn*Dn+Knn),Hinv,Hinv*(Cnn*Enn+C21*Enc+M21*Bnc);...
            Hconst*Hinv*(Cnn*Dn+Knn)+C12*Dn+K12+M12*An,Hconst*Hinv,M12*Bnn+C11*Enc+C12*Enn+...
            Hconst*Hinv*(M21*Bnc+C21*Enc+Cnn*Enn);...
            zeros(1,4),1];
    elseif i == SegNum+1
        Fconst = C22*Dn+K22;
        Hconst = C22*Dn+K22;
        P(:,:,i-1) = [eye(2) zeros(2,3);M22*An eye(2) M22*Bnc;zeros(1,4),1];
        F(:,:,i-1) = [Finv*(Cnn*Dn+K11),Finv,Finv*(C12*EnL+C12*Enc+M12*Bnc);...
            Fconst*Finv*(Cnn*Dn+Knn)+C21*Dn+K21+M21*An,Fconst*Finv,C22*Enc+C21*EnL+M21*BnL+...
            Fconst*Finv*(Cnn*EnL+C12*Enc+M12*Bnc);...
            zeros(1,4),1];
        J(:,:,i-1) = [eye(2) zeros(2,3);-M22*An eye(2) -M22*Bnc;zeros(1,4),1];
        H(:,:,i-1) = [Hinv*(C22*Dn+K22),Hinv,Hinv*(C22*Enn+C21*Enc+M21*Bnc);...
            Hconst*Hinv*(C22*Dn+K22)+C12*Dn+K12+M12*An,Hconst*Hinv,M12*Bnn+Cnn*Enc+C12*Enn+...
            Hconst*Hinv*(M21*Bnc+C21*Enc+C22*Enn);...
            zeros(1,4),1];
    else
        Fconst = Cnn*Dn+Knn;
        Hconst = Cnn*Dn+Knn;
        P(:,:,i-1) = [eye(2) zeros(2,3);Mnn*An eye(2) Mnn*Bnc;zeros(1,4),1];
        F(:,:,i-1) = [Finv*(Cnn*Dn+Knn),Finv,Finv*(C12*EnL+C12*Enc+M12*Bnc);...
            Fconst*Finv*(Cnn*Dn+Knn)+C21*Dn+K21+M21*An,Fconst*Finv,Cnn*Enc+C21*EnL+M21*BnL+...
            Fconst*Finv*(Cnn*EnL+C12*Enc+M12*Bnc);...
            zeros(1,4),1];
        J(:,:,i-1) = [eye(2) zeros(2,3);-Mnn*An eye(2) -Mnn*Bnc;zeros(1,4),1];
        H(:,:,i-1) = [Hinv*(Cnn*Dn+Knn),Hinv,Hinv*(Cnn*Enn+C21*Enc+M21*Bnc);...
            Hconst*Hinv*(Cnn*Dn+Knn)+C12*Dn+K12+M12*An,Hconst*Hinv,M12*Bnn+Cnn*Enc+C12*Enn+...
            Hconst*Hinv*(M21*Bnc+C21*Enc+Cnn*Enn);...
            zeros(1,4),1];
    end
end

%Use computed P, F, H, and J to create Q and T to generate state and output
%matrix
%Note: K.. C.. and M.. are recalculated because I am lazy. Fix later
for j = node-1:node+1
    K11 = E*I/L(j-1)^3*[12 6;6 4];
    K12 = E*I/L(j-1)^3*[-12 6;-6 2];
    K21 = E*I/L(j-1)^3*[-12 -6;6 2];
    K22 = E*I/L(j-1)^3*[12 -6;-6 4];
    Knn = K11+K22;

    C11 = etaE*I/L(j-1)^3*[12 6;6 4];
    C12 = etaE*I/L(j-1)^3*[-12 6;-6 2];
    C21 = etaE*I/L(j-1)^3*[-12 -6;6 2];
    C22 = etaE*I/L(j-1)^3*[12 -6;-6 4];
    Cnn = C11+C22;

    M11 = rho*CA/(420*L(j-1))*[156 22;22 4];
    M12 = rho*CA/(420*L(j-1))*[54 -13;13 -3];
    M21 = rho*CA/(420*L(j-1))*[54 13;-13 -3];
    M22 = rho*CA/(420*L(j-1))*[156 -22;-22 4];
    Mnn = M11+M22;
    
    Q = eye(5);
    for i = 1:j-1
        Q = P(:,:,i)*F(:,:,i)*Q;
    end

    T = eye(5);
    for i = fliplr(j+2:SegNum)
        T = H(:,:,i)*J(:,:,i)*T;
    end
    T = J(:,:,node+1)*T;
    
    CXn = Mnn*An-(T(3:4,1:2)*T(1:2,1:2)\F(1:2,3:4,node+1)-F(3:4,3:4,node+1))\...
        (F(3:4,1:2,node+1)-T(3:4,1:2)*T(1:2,1:2)\F(1:2,1:2,node+1))+...
        (Q(3:4,3:4)*Q(1:2,3:4)\H(1:2,3:4,node)-H(3:4,3:4,node))\(H(3:4,1:2,node)-Q(3:4,3:4)*Q(1:2,3:4)\H(1:2,1:2,node));

    CBn = (T(3:4,1:2)*T(1:2,1:2)\F(1:2,3:4,node+1)-F(3:4,3:4,node+1))\Fconst*Finv*M21...
        -(Q(3:4,3:4)*Q(1:2,3:4)\H(1:2,3:4)-H(3:4,3:4))\Hconst*Hinv*M12-Mnn;

    CEn = (T(3:4,1:2)*T(1:2,1:2)\F(1:2,3:4,node+1)-F(3:4,3:4,node+1))\(Fconst*Finv*Cnn+C21-T(3:4,1:2)*T(1:2,1:2)\(Finv*Cnn))...
        -(Q(3:4,3:4)*Q(1:2,3:4)\H(1:2,3:4)-H(3:4,3:4))\(Hconst*Hinv*Cnn-C12+Q(3:4,3:4)*Q(1:2,3:4)\(Hinv*Cnn));

    CEnp1 = (T(3:4,1:2)*T(1:2,1:2)\F(1:2,3:4,node+1)-F(3:4,3:4,node+1))\...
        (Fconst*Finv*C12+Cnn-T(3:4,1:2)*T(1:2,1:2)\(Finv*C12)+T(3:4,1:2)*T(1:2,1:2)\(Hinv*C21)+Mnn*An*Hinv*C21-(Hconst*Hinv*C21+Cnn));

    CEnm1 = -(Q(3:4,3:4)*Q(1:2,3:4)\H(1:2,3:4)-H(3:4,3:4,node))\...
        (Hconst*Hinv*C21+Cnn-Q(3:4,3:4)*Q(1:2,3:4)\(Hinv*C21)-Mnn*An*Finv*C12-Fconst*Finv*C12-Cnn+Q(3:4,3:4)*Q(1:2,3:4)\(Finv*C12));

    CBnp1 = (T(3:4,1:2)*T(1:2,1:2)\F(1:2,3:4,node+1)-F(3:4,3:4,node+1))\...
        (M12-T(3:4,1:2)*T(1:2,1:2)\(Finv*M12)+T(3:4,1:2)*T(1:2,1:2)\(Hinv*M21)+Mnn*An*Hinv*M21-M21+Mnn);

    CBnm1 = -(Q(3:4,3:4)*Q(1:2,3:4)\H(1:2,3:4)-H(3:4,3:4,node))\...
        (M21-Q(3:4,3:4)*Q(1:2,3:4)\(Hinv*M21)+Q(3:4,3:4)*Q(1:2,3:4)\(Finv*M12)-Mnn*An*Finv*M12-M12-Mnn);
    
    A((j-node+1)*6+1:(j-node+2)*6,(j-node+1)*6+1:(j-node+1)*6+18) = [CXn\(-3*CEnm1/dt+5*CBnm1/dt^2),CXn\(3*CEnm1/(2*dt)-4*CBnm1/dt^2),CXn\(-CEnm1/(3*dt)+CBnm1/dt^2),...
        CXn\(-3*CEn/dt+5*CBn/dt^2),CXn\(3*CEn/(2*dt)-4*CBn/dt^2),CXn\(-CEn/(3*dt)+CBn/dt^2),...
        CXn\(-3*CEnp1/dt+5*CBnp1/dt^2),CXn\(3*CEnp1/(2*dt)-4*CBnp1/dt^2),CXn\(-CEnp1/(3*dt)+CBnp1/dt^2);
        zeros(4,6),[eye(2),zeros(2,4);zeros(2),eye(2),zeros(2)],zeros(4,6)];
    B((j-node+1)*6+1:(j-node+2)*6,(j-node+1)*2+1:(j-node+1)*2+2) = [CXn\eye(2);zeros(4,2)];
end

A = A(:,7:24);
B = B(:,3:4);