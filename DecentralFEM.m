function [A,B] = DecentralFEM(dt,E,etaE,I,rho,CA,L,x_0,dx_0,ddx,x2,ddx2,x3,int,node)

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
P = zeros(5,5,length(x_0)-1);
F = zeros(5,5,length(x_0)-1);
H = zeros(5,5,length(x_0)-1);
J = zeros(5,5,length(x_0)-1);

K11 = E*I/L^3*[12 6;6 4];
K12 = E*I/L^3*[-12 6;-6 2];
K21 = E*I/L^3*[-12 -6;6 2];
K22 = E*I/L^3*[12 -6;-6 4];

C11 = etaE*I/L^3*[12 6;6 4];
C12 = etaE*I/L^3*[-12 6;-6 2];
C21 = etaE*I/L^3*[-12 -6;6 2];
C22 = etaE*I/L^3*[12 -6;-6 4];

M11 = rho*CA/(420*L)*[156 22;22 4];
M12 = rho*CA/(420*L)*[54 -13;13 -3];
M21 = rho*CA/(420*L)*[54 13;-13 -3];
M22 = rho*CA/(420*L)*[156 -22;-22 4];
Mnn = M11+M22;

[An,Bn,Dn,En] = Houbolt(dt,x_0,dx_0,ddx,x2,ddx2,x3,int);
Finv = pinv(-C12*Dn-K12-M12*An);
Fconst = C22*Dn+K22;
Hinv = pinv(-C21*Dn-K21-M21*An);
Hconst = C22*Dn+K22;

for i = 2:SegNum+1
    BnL = Bn(:,i-1);
    EnL = En(:,i-1);
    Bnc = Bn(:,i);
    Enc = En(:,i);
    Bnn = Bn(:,i+1);
    Enn = En(:,i+1);
    P(:,:,i-1) = [eye(2) zeros(2,3);Mnn*An eye(2) Mnn*Bnc;zeros(1,4),1];
    F(:,:,i-1) = [Finv*(C11*Dn+K11),Finv,Finv*(C12*EnL+C12*Enc+M12*Bnc);...
        Fconst*Finv*(C11*Dn+K11)+C21*Dn+K21+M21*An,Fconst*Finv,C22*Enc+C21*EnL+M21*BnL+...
        Fconst*Finv*(C11*EnL+C12*Enc+M12*Bnc);...
        zeros(1,4),1];
    J(:,:,i-1) = [eye(2) zeros(2,3);-Mnn*An eye(2) -Mnn*Bnc;zeros(1,4),1];
    H(:,:,i-1) = [Hinv*(C22*Dn+K22),Hinv,Hinv*(C22*Enn+C21*Enc+M21*Bnc);...
        Hconst*Hinv*(C22*Dn+K22)+C12*Dn+K12+M12*An,Hconst*Hinv,M12*Bnn+C11*Enc+C12*Enn+...
        Hconst*Hinv*(M21*Bnc+C21*Enc+C22*Enn);...
        zeros(1,4),1];
end

Q = eye(5);
for i = 1:node-1
        Q = P(:,:,i)*F(:,:,i)*Q;
end

G = eye(5);
T = eye(5);
for i = fliplr(node+2:SegNum)
    T = J(:,:,i)*H(:,:,i)*T;
end
R = eye(5);

CXn = Mnn*An-(T(3:4,1:2)*(T(1:2,1:2)\(R(1:2,1:2)*F(1:2,3:4,node+1)+R(1:2,3:4)*F(3:4,3:4,node+1)))-...
    (R(3:4,1:2)*F(1:2,3:4,node+1)+R(3:4,3:4)*F(3:4,3:4,node+1)))\...
    (-T(3:4,1:2)*(T(1:2,1:2)\(R(1:2,1:2)*F(1:2,1:2,node+1)+R(1:2,3:4)*F(3:4,1:2,node+1)))+...
    R(3:4,1:2)*F(1:2,1:2,node+1)+R(3:4,3:4)*F(3:4,1:2,node+1))+...
    (Q(3:4,3:4)*(Q(1:2,3:4)\(G(1:2,1:2)*H(1:2,3:4,node)+G(1:2,3:4)*H(3:4,3:4,node)))-...
    (G(3:4,1:2)*H(1:2,3:4,node)+G(3:4,3:4)*H(3:4,3:4,node)))\...
    (G(3:4,1:2)*H(1:2,1:2,node)+G(3:4,3:4)*H(3:4,1:2,node)-Q(3:4,3:4)*...
    (Q(1:2,3:4)\(G(1:2,1:2)*H(1:2,1:2,node)+G(1:2,3:4)*H(3:4,1:2,node))));

CBn = (T(3:4,1:2)*(T(1:2,1:2)\(R(1:2,1:2)*F(1:2,3:4,node+1)+R(1:2,3:4)*F(3:4,1:2,node+1)))-...
    (R(3:4,1:2)*F(1:2,3:4,node+1)+R(3:4,3:4)*F(3:4,3:4,node+1)))\...
    ((R(3:4,3:4)*Fconst*Finv-T(3:4,1:2)*(T(1:2,1:2)\(R(1:2,1:2)*Finv+R(1:2,3:4)*Fconst*Finv))+...
    R(1:2,1:2)*Finv)*M21)-Mnn-(Q(3:4,3:4)*(Q(1:2,3:4)\(G(1:2,1:2)*H(1:2,3:4,node)+G(1:2,3:4)*H(3:4,3:4,node)))-...
    (G(3:4,1:2)*H(1:2,3:4,node)+G(3:4,3:4)*H(3:4,3:4,node)))\...
    (-Q(3:4,3:4)*(Q(1:2,3:4)\(G(1:2,1:2)*Hinv+G(1:2,3:4)*Hconst*Hinv))+...
    G(3:4,3:4)*Hconst*Hinv+G(3:4,1:2)*Hinv)*M12;

CEn = (T(3:4,1:2)*(T(1:2,1:2)\(R(1:2,1:2)*F(1:2,3:4,node+1)+R(1:2,3:4)*F(3:4,1:2,node+1)))-...
    (R(3:4,1:2)*F(1:2,3:4,node+1)+R(3:4,3:4)*F(3:4,3:4,node+1)))\...
    (R(3:4,1:2)*Finv*C11+R(3:4,3:4)*(Fconst*Finv*C11+C21)-T(3:4,1:2)*(T(1:2,1:2)\...
    (R(1:2,1:2)*Finv*C11+R(1:2,3:4)*(Fconst*Finv*C11+C21))))...
    -(Q(3:4,3:4)*(Q(1:2,3:4)\(G(1:2,1:2)*H(1:2,3:4,node)+G(1:2,3:4)*H(3:4,3:4,node)))-...
    (G(3:4,1:2)*H(1:2,3:4,node)+G(3:4,3:4)*H(3:4,3:4,node)))\...
    (-Q(3:4,3:4)*(Q(1:2,3:4)\(G(1:2,1:2)*Hinv*C22+G(1:2,3:4)*(Hconst*Hinv*C22+C12)))+...
    G(3:4,1:2)*Hinv*C22+G(3:4,3:4)*(Hconst*Hinv*C22+C12));

A = [CXn\(-3*CEn/dt+5*CBn/dt^2),CXn\(3*CEn/(2*dt)-4*CBn/dt^2),CXn\(-CEn/(3*dt)+CBn/dt^2);
    eye(2),zeros(2,4);zeros(2),eye(2),zeros(2)];

B = [pinv(CXn);zeros(4,2)];