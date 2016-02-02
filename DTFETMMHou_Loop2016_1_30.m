function [x,e] = DTFETMMHou_Loop2015_1_30(dt,E,etaE,I,rho,CA,L,xL_0,dxL_0,ddxL_0,xR_0,dxR_0,ddxR_0,fn,Tf,t_imp)
t = cputime;
T = 0:dt:Tf;
xL = zeros(4,length(xL_0)+1,length(T));
dxL = zeros(2,length(xL_0)+1,length(T));
ddxL = zeros(2,length(xL_0)+1,length(T));
xL(:,:,1) = [zeros(4,1),[x_0;zeros(2,length(x_0))]];
dxL(:,:,1) = [zeros(2,1),dx_0];
ddxL(:,:,1) = [zeros(2,1),ddx_0];

for i = 1:length(T)
    if T(i) < t_imp
        if i == 1
            xL(:,:,i+1) = DTTMMHou2015_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,1),dxL(:,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),i,fn);
        elseif i == 2
            xL(:,:,i+1) = DTTMMHou2015_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-1),i,fn);
        else
            xL(:,:,i+1) = DTTMMHou2015_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),x(1:2,2:end,i-2),i,fn);
        end
    else
        if i == 1
            xL(:,:,i+1) = DTTMMHou2015_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,1),dxL(:,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),i,zeros(2,20));
        elseif i == 2
            xL(:,:,i+1) = DTTMMHou2015_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-1),i,zeros(2,20));
        else
            xL(:,:,i+1) = DTTMMHou2015_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-2),i,zeros(2,20));
        end
    end
    if i == 1
        [An,Bn,Dn,En] = Houbolt(dt,xL(1:2,2:end,1),dxL(:,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),i);
    elseif i == 2
        [An,Bn,Dn,En] = Houbolt(dt,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i),i);
    else
        [An,Bn,Dn,En] = Houbolt(dt,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-2),i);
    end
    dx(:,:,i+1) = Dn*x(1:2,:,i+1)+[[0;0],En];
    ddx(:,:,i+1) = An*x(1:2,:,i+1)+[[0;0],Bn];
end
e = cputime-t;