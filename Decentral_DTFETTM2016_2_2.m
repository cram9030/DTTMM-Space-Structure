%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = Decentral_DTFETTM2016_2_2(dt,E,etaE,I,rho,CA,L,SegNum,node)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decentral_DTFETTM2016_1_31: Create decentralized model of space structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is called to create an decentralize model of a two node
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

K11 = zeros(2,2,length(x_0)-1);
K12 = zeros(2,2,length(x_0)-1);
K21 = zeros(2,2,length(x_0)-1);
K22 = zeros(2,2,length(x_0)-1);
Knn = zeros(2,2,length(x_0)-1);

C11 = zeros(2,2,length(x_0)-1);
C12 = zeros(2,2,length(x_0)-1);
C21 = zeros(2,2,length(x_0)-1);
C22 = zeros(2,2,length(x_0)-1);
Cnn = zeros(2,2,length(x_0)-1);

M11 = zeros(2,2,length(x_0)-1);
M12 = zeros(2,2,length(x_0)-1);
M21 = zeros(2,2,length(x_0)-1);
M22 = zeros(2,2,length(x_0)-1);
Mnn = zeros(2,2,length(x_0)-1);

Fconst = zeros(2,2,length(x_0)-1);
Finv = zeros(2,2,length(x_0)-1);
Hconst = zeros(2,2,length(x_0)-1);
Hinv = zeros(2,2,length(x_0)-1);

%Generate P, F, H, and J
for i = 2:SegNum+1
    K11(:,:,i-1) = E*I/L(i-1)^3*[12 6;6 4];
    K12(:,:,i-1) = E*I/L(i-1)^3*[-12 6;-6 2];
    K21(:,:,i-1) = E*I/L(i-1)^3*[-12 -6;6 2];
    K22(:,:,i-1) = E*I/L(i-1)^3*[12 -6;-6 4];
    Knn(:,:,i-1) = K11(:,:,i-1)+K22(:,:,i-1);

    C11(:,:,i-1) = etaE*I/L(i-1)^3*[12 6;6 4];
    C12(:,:,i-1) = etaE*I/L(i-1)^3*[-12 6;-6 2];
    C21(:,:,i-1) = etaE*I/L(i-1)^3*[-12 -6;6 2];
    C22(:,:,i-1) = etaE*I/L(i-1)^3*[12 -6;-6 4];
    Cnn(:,:,i-1) = C11(:,:,i-1)+C22(:,:,i-1);

    M11(:,:,i-1) = rho*CA/(420*L(i-1))*[156 22;22 4];
    M12(:,:,i-1) = rho*CA/(420*L(i-1))*[54 -13;13 -3];
    M21(:,:,i-1) = rho*CA/(420*L(i-1))*[54 13;-13 -3];
    M22(:,:,i-1) = rho*CA/(420*L(i-1))*[156 -22;-22 4];
    Mnn(:,:,i-1) = M11(:,:,i-1)+M22(:,:,i-1);

    [An,Bn,Dn,En] = Houbolt(dt,x_0,dx_0,ddx,x2,ddx2,x3,int);
    Finv(:,:,i-1) = (-C12(:,:,i-1)*Dn-K12(:,:,i-1)-M12(:,:,i-1)*An)\eye(2);
    Hinv(:,:,i-1) = (-C21(:,:,i-1)*Dn-K21(:,:,i-1)-M21(:,:,i-1)*An)\eye(2);

    BnL = Bn(:,i-1);
    EnL = En(:,i-1);
    Bnc = Bn(:,i);
    Enc = En(:,i);
    Bnn = Bn(:,i+1);
    Enn = En(:,i+1);
    if i == 2
        Fconst(:,:,i-1) = C11(:,:,i-1)*Dn+K11(:,:,i-1);
        Hconst(:,:,i-1) = C11(:,:,i-1)*Dn+K11(:,:,i-1);
        P(:,:,i-1) = [eye(2) zeros(2,3);M11(:,:,i-1)*An eye(2) M11(:,:,i-1)*Bnc;zeros(1,4),1];
        F(:,:,i-1) = [Finv*(C11(:,:,i-1)*Dn+K11(:,:,i-1)),Finv,Finv*(C12(:,:,i-1)*EnL+C12(:,:,i-1)*Enc+M12(:,:,i-1)*Bnc);...
            Fconst(:,:,i-1)*Finv(:,:,i-1)*(C11(:,:,i-1)*Dn+K11(:,:,i-1))+C21(:,:,i-1)*Dn+K21(:,:,i-1)+M21(:,:,i-1)*An,Fconst(:,:,i-1)*Finv(:,:,i-1),Cnn(:,:,i-1)*Enc+C21(:,:,i-1)*EnL+M21(:,:,i-1)*BnL+...
            Fconst(:,:,i-1)*Finv(:,:,i-1)*(C11(:,:,i-1)*EnL+C12(:,:,i-1)*Enc+M12(:,:,i-1)*Bnc);...
            zeros(1,4),1];
        J(:,:,i-1) = [eye(2) zeros(2,3);-M11(:,:,i-1)*An eye(2) -M11(:,:,i-1)*Bnc;zeros(1,4),1];
        H(:,:,i-1) = [Hinv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+Knn(:,:,i-1)),Hinv(:,:,i-1),Hinv(:,:,i-1)*(Cnn(:,:,i-1)*Enn+C21(:,:,i-1)*Enc+M21(:,:,i-1)*Bnc);...
            Hconst(:,:,i-1)*Hinv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+Knn(:,:,i-1))+C12(:,:,i-1)*Dn+K12(:,:,i-1)+M12(:,:,i-1)*An,Hconst(:,:,i-1)*Hinv(:,:,i-1),M12(:,:,i-1)*Bnn+C11(:,:,i-1)*Enc+C12(:,:,i-1)*Enn+...
            Hconst(:,:,i-1)*Hinv(:,:,i-1)*(M21(:,:,i-1)*Bnc+C21(:,:,i-1)*Enc+Cnn(:,:,i-1)*Enn);...
            zeros(1,4),1];
    elseif i == SegNum+1
        Fconst(:,:,i-1) = C22(:,:,i-1)*Dn+K22(:,:,i-1);
        Hconst(:,:,i-1) = C22(:,:,i-1)*Dn+K22(:,:,i-1);
        P(:,:,i-1) = [eye(2) zeros(2,3);M22(:,:,i-1)*An eye(2) M22(:,:,i-1)*Bnc;zeros(1,4),1];
        F(:,:,i-1) = [Finv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+K11(:,:,i-1)),Finv(:,:,i-1),Finv(:,:,i-1)*(C12(:,:,i-1)*EnL+C12(:,:,i-1)*Enc+M12(:,:,i-1)*Bnc);...
            Fconst(:,:,i-1)*Finv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+Knn(:,:,i-1))+C21(:,:,i-1)*Dn+K21(:,:,i-1)+M21(:,:,i-1)*An,Fconst*Finv,C22(:,:,i-1)*Enc+C21(:,:,i-1)*EnL+M21(:,:,i-1)*BnL+...
            Fconst(:,:,i-1)*Finv(:,:,i-1)*(Cnn(:,:,i-1)*EnL+C12(:,:,i-1)*Enc+M12(:,:,i-1)*Bnc);...
            zeros(1,4),1];
        J(:,:,i-1) = [eye(2) zeros(2,3);-M22(:,:,i-1)*An eye(2) -M22(:,:,i-1)*Bnc;zeros(1,4),1];
        H(:,:,i-1) = [Hinv(:,:,i-1)*(C22(:,:,i-1)*Dn+K22(:,:,i-1)),Hinv(:,:,i-1),Hinv(:,:,i-1)*(C22(:,:,i-1)*Enn+C21(:,:,i-1)*Enc+M21(:,:,i-1)*Bnc);...
            Hconst(:,:,i-1)*Hinv(:,:,i-1)*(C22(:,:,i-1)*Dn+K22(:,:,i-1))+C12(:,:,i-1)*Dn+K12(:,:,i-1)+M12(:,:,i-1)*An,Hconst(:,:,i-1)*Hinv(:,:,i-1),M12(:,:,i-1)*Bnn+Cnn(:,:,i-1)*Enc+C12(:,:,i-1)*Enn+...
            Hconst(:,:,i-1)*Hinv(:,:,i-1)*(M21(:,:,i-1)*Bnc+C21(:,:,i-1)*Enc+C22(:,:,i-1)*Enn);...
            zeros(1,4),1];
    else
        Fconst(:,:,i-1) = Cnn(:,:,i-1)*Dn+Knn(:,:,i-1);
        Hconst(:,:,i-1) = Cnn(:,:,i-1)*Dn+Knn(:,:,i-1);
        P(:,:,i-1) = [eye(2) zeros(2,3);Mnn(:,:,i-1)*An eye(2) Mnn(:,:,i-1)*Bnc;zeros(1,4),1];
        F(:,:,i-1) = [Finv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+Knn(:,:,i-1)),Finv(:,:,i-1),Finv(:,:,i-1)*(C12(:,:,i-1)*EnL+C12(:,:,i-1)*Enc+M12(:,:,i-1)*Bnc);...
            Fconst(:,:,i-1)*Finv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+Knn(:,:,i-1))+C21(:,:,i-1)*Dn+K21(:,:,i-1)+M21(:,:,i-1)*An,Fconst(:,:,i-1)*Finv(:,:,i-1),Cnn(:,:,i-1)*Enc+C21(:,:,i-1)*EnL+M21(:,:,i-1)*BnL+...
            Fconst(:,:,i-1)*Finv(:,:,i-1)*(Cnn(:,:,i-1)*EnL+C12(:,:,i-1)*Enc+M12(:,:,i-1)*Bnc);...
            zeros(1,4),1];
        J(:,:,i-1) = [eye(2) zeros(2,3);-Mnn(:,:,i-1)*An eye(2) -Mnn(:,:,i-1)*Bnc;zeros(1,4),1];
        H(:,:,i-1) = [Hinv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+Knn(:,:,i-1)),Hinv(:,:,i-1),Hinv(:,:,i-1)*(Cnn(:,:,i-1)*Enn+C21(:,:,i-1)*Enc+M21(:,:,i-1)*Bnc);...
            Hconst(:,:,i-1)*Hinv(:,:,i-1)*(Cnn(:,:,i-1)*Dn+Knn(:,:,i-1))+C12(:,:,i-1)*Dn+K12(:,:,i-1)+M12(:,:,i-1)*An,Hconst(:,:,i-1)*Hinv(:,:,i-1),M12(:,:,i-1)*Bnn+Cnn(:,:,i-1)*Enc+C12(:,:,i-1)*Enn+...
            Hconst(:,:,i-1)*Hinv(:,:,i-1)*(M21(:,:,i-1)*Bnc+C21(:,:,i-1)*Enc+Cnn(:,:,i-1)*Enn);...
            zeros(1,4),1];
    end
end

%Use computed P, F, H, and J to create Q and T to generate state and output
%matrix

Q = eye(5);
for i = 1:node-1
    Q = P(:,:,i)*F(:,:,i)*Q;
end

T = eye(5);
for i = fliplr(node+2:SegNum)
    T = H(:,:,i)*J(:,:,i)*T;
end
T = J(:,:,node+1)*T;

CnXn = Mnn(:,:,node)*An-(T(1:2,1:2)*T(1:2,3:4)\F(3:4,3:4,node+1)-F(1:2,3:4,node+1))\...
    (F(1:2,1:2,node+1)-T(1:2,1:2)*T(1:2,3:4)\F(3:4,1:2,node+1))-H(1:2,3:4,node)\H(1:2,1:2,node);

CnXnm1 = H(1:2,1:2,node)\eye(2);

CnBn = -(T(1:2,1:2)*T(1:2,3:4)\F(3:4,3:4,node+1)-F(1:2,3:4,node+1))\...
    (T(1:2,1:2)*T(1:2,3:4)\(Fconst(:,:,node+1)*Finv(:,:,node+1))*M21(:,:,node))+Mnn(:,:,node);

CnBnm1 = H(1:2,3:4,node)\Hinv(:,:,node)*M21(:,:,node-1);

CnEnm1 = H(1:2,3:4,node)\Hinv(:,:,node)*C21(:,:,node-1);

CnEn = (T(1:2,1:2)*T(1:2,3:4)\F(3:4,3:4,node+1)-F(1:2,3:4,node+1))\...
    (Finv(:,:,node+1)*Cnn(:,:,node)-T(1:2,1:2)*T(1:2,3:4)\((Fconst(:,:,node-1)*Finv(:,:,node-1))+...
    Cnn(:,:,node)+C21(:,:,node)))+H(1:2,3:4,node-1)*Hconst(:,:,node)*Cnn(:,:,node);

CnBnp1 = (T(1:2,1:2)*T(1:2,3:4)\F(3:4,3:4,node+1)-F(1:2,3:4,node+1))\...
    (Finv(:,:,node+1)*M12(:,:,node+1)-T(1:2,1:2)*T(1:2,3:4)\M12(:,:,node+1));

CnEnp1 = (T(1:2,1:2)*T(1:2,3:4)\F(3:4,3:4,node+1)-F(1:2,3:4,node+1))\...
    (Finv(:,:,node)*C12(:,:,node+1)-T(1:2,1:2)*T(1:2,3:4)\(Fconst(:,:,node+1)*Finv(:,:,node+1)*C12(:,:,node+1)+Cnn(:,:,node+1)));

Cnode = (T(1:2,1:2)*T(1:2,3:4)\F(3:4,3:4,node+1)-F(1:2,3:4,node+1))\...
    (T(1:2,1:2)*T(1:2,3:4)\T(3:4,5)-T(1:2,5));

Cnm1Xnm1 = Mnn(:,:,node-1)*An+(Q(1:2,3:4)*Q(3:4,3:4)\H(3:4,3:4,node-1)-H(1:2,3:4,node-1))\...
    (H(1:2,1:2,node-1)-Q(1:2,3:4)*Q(3:4,3:4)\H(3:4,1:2,node-1))-H(3:4,3:4,node)*H(1:2,1:2,node)\eye(2);

Cnm1Xn = H(3:4,1:2,node)-H(3:4,3:4,node)*H(1:2,3:4,node)\H(1:2,1:2,node);

Cnm1Bn = Hconst(:,:,node)*Hinv(:,:,node)*M12(:,:,node);

Cnm1Bnm1 = M21(:,:,node-1)-H(3:4,3:4,node)*H(1:2,3:4,node)\Hinv(:,:,node)*M21(:,:,node-1)+...
    (Q(1:2,3:4)*Q(3:4,3:4)\H(3:4,3:4,node-1)-H(1:2,3:4,node-1))\(Q(1:2,3:4)*Q(3:4,3:4)\...
    Hconst(:,:,node-1)*Hinv(:,:,node-1)*M12(:,:,node-1))-Mnn(:,:,node-1);

Cnm1Bnm2 = -(Hinv(:,:,node-1)*M21(:,:,node-2)-Q(1:2,3:4)*Q(3:4,3:4)\H(3:4,1:2,node-1));

Cnm1Enm1 = Hconst(:,:,node)*Hinv(:,:,node)*C21(:,:,node-1)+Cnn(:,:,node-1)-H(3:4,3:4,node)*H(1:2,3:4,node)\Hinv(:,:,node)*C21(:,:,node-1)-...
    (Q(1:2,3:4)*Q(3:4,3:4)\H(3:4,3:4,node-1)-H(1:2,3:4,node-1))\(Hinv(:,:,node-1)*Cnn(:,:,node-1)-Q(1:2,3:4)*Q(3:4,3:4)\...
    (Hconst(:,:,node-1)*Hinv(:,:,node-1)*Cnn(:,:,node-1)+C12(:,:,node-1)));

Cnm1En = Hconst(:,:,node)*Hinv(:,:,node)*Cnn(:,:,node)+C12(:,:,node)-H(3:4,3:4,node)*H(1:2,3:4,node)\Hinv(:,:,node)*Cnn(:,:,node);

Cnm1Enm2 = (Q(1:2,3:4)*Q(3:4,3:4)\H(3:4,3:4,node-1)-H(1:2,3:4,node-1))\...
    (Hinv(:,:,node-1)*Cnn(:,:,node-2)-Q(1:2,3:4)*Q(3:4,3:4)\...
    (Hconst(:,:,node-1)*Hinv(:,:,node-1)*C21(:,:,node-2)+Cnn(:,:,node-2)));

Cnodem1 = -(Q(1:2,3:4)*Q(3:4,3:4)\H(3:4,3:4,node-1)-H(1:2,3:4,node-1))\...
    (Q(1:2,3:4)*Q(3:4,3:4)\Q(3:4,5)-Q(1:2,5));

StateMatrix = [Cnm1Xnm1,Cnm1Xn;CnXnm1,CnXn];
AccMatrix = [Cnm1Bnm2,Cnm1Bn1,Cnm1Bn,zeros(2);zeros(2),CnBnm1,CnBn,CnBnp1];
VelMatrix = [Cnm1Enm2,Cnm1En1,Cnm1En,zeros(2);zeros(2),CnEnm1,CnBn,CnEnp1];

A = 