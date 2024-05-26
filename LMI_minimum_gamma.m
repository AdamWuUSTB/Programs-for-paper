clear;clc;
global A;global B;global C; global D;global n;global m;global p;global q;
A=[-2.98 0.4 0 -0.034;-0.99 -0.21 0.035 -0.0011;0 0 0 1;0.39 -5.55 0 -1.89];
B=[-0.032;0;0;-1.6];D=[0.5;0;0.6;0];%飞行器参数取0.4
C=[0 0 1 0];
%%%%%%构造误差系统
A0=[zeros(p,p) C;zeros(n,p) A];
B0=[zeros(p,m);B];
C0=[eye(p,p) zeros(p,n)];
D0=[-eye(p,p) zeros(p,q);zeros(n,p)  D];
%%%%%%求解线性矩阵不等式
gamma=1.5;q_e=1.5; q_x=0.1; r=0.1;
M=[sqrt(q_e*eye(p,p)) zeros(p,n);zeros(n,p) sqrt(q_x*eye(n,n));zeros(m,p) zeros(m,n)];
N=[zeros(p,m);zeros(n,m);sqrt(r*eye(m,m))];
%%%初始化LMI
setlmis([]);
%%%定义决策变量
U=lmivar(1,[m,1]);
P=lmivar(1,[p+n,1]);
H=lmivar(2,[m,n+p]);
L=lmivar(2,[m,n+p]);
KK=lmivar(2,[1,1]);
%%%LMI的描述
% 1st LMI
lmiterm([1 1 1 P],A0,1,'s');
lmiterm([1 1 1 L],B0,1,'s');
lmiterm([1 1 3 0],D0);
lmiterm([1 2 1 U],-1,B0');
lmiterm([1 2 1 H],1,1);
lmiterm([1 2 1 L],1,1);
lmiterm([1 2 2 U],-2,1);
lmiterm ([1 3 3 KK],-1,1);
lmiterm([1 4 1 P],M,1);
lmiterm([1 4 1 L],N,1);
lmiterm([1 4 2 U],N,-1);
lmiterm([1 4 4 0],-1);
%2rd LMI
    for i=1:m
        lmiterm([-2 1 1 P],1,1);
        lmiterm([-2 2 1 H(i,:)],1,1);
        lmiterm([-2 2 2 0],25);
    end
% % 3rd-4th LMI
lmiterm([-4 1 1 P],1,1);
lmiterm([-5 1 1 U],1,1);
%%%LMI的求解
% lmis=getlmis;
% [tmin,xfeas]=feasp(lmis);
% U=dec2mat(lmis,xfeas,U);
% P=dec2mat(lmis,xfeas,P);
% H=dec2mat(lmis,xfeas,H);
% L=dec2mat(lmis,xfeas,L);
%控制增益阵的求解;

% K=L*inv(P)
% [k11,k12]=size(K);
% K_e=K(:,1:q)
% K_x=K(:,q+1:k12)

%%%%%%%%%%% 4 %%%%%%%%%%
lmiterm([-6 1 1 KK],1,1);
lmisys=getlmis;

nnn=decnbr(lmisys);
c=zeros(1,nnn);
for j=1:nnn
    [KKj]=defcx(lmisys,j,KK);
    c(j)=KKj;
end

[copt,xopt]=mincx(lmisys,c);

%L11=L
%U11=U
%X11=X
% lmiinfo(lmisys);
U=dec2mat(lmisys,xopt,U);
P=dec2mat(lmisys,xopt,P);
H=dec2mat(lmisys,xopt,H);
L=dec2mat(lmisys,xopt,L);
KK=dec2mat(lmisys,xopt,KK);

K=L*inv(P);
gamma=KK^(1/2)



