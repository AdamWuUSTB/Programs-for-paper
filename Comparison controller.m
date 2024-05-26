function SPC
clear;clc;
global A;global Bd;global B;global C;global tf;global s;global D1;global B1;global C1;
global n;global m;global p;global q;global h;global llr;global u0;
global M1; global M2; global Ec;
llr=1.00;llr1=0.50;%设定预见长度
A=[-2.98 0.4 0 -0.034;-0.99 -0.21 0.035 -0.0011;0 0 0 1;0.39 -5.55 0 -1.89];
B=[-0.032;0;0;-1.6];Bd=[0.5;0;0.6;0];%飞行器参数取0.4
C=[0 0 1 0];
M1=[0 0 0 1;1 0 0 0;0 1 0 0];M2=[C;M1]; 
% A=[0 1;0 0];
% B=[0;1];Bd=[1;1];
% C=[1 0];
% M1=[0 1];M2=[C;M1]; 
u0=5;
[n,n]=size(A); [n,m]=size(B); 
[p,n]=size(C); [n,q]=size(Bd);
t0=0; tf=50; h=0.001; s=t0:h:tf; [k0,k]=size(s);ss=t0:h:tf+llr+llr;[k0,k1]=size(ss);
[U,P,L,H,K_e,K_x] = LMI(k);
[F_x,F_e,F_ksai,Ec]=LMI1(k);

for i=1:k  
    ww(i)=disturbance_2(s(i));
end
for i=1:k
    lr(i)=ref_signal_2(s(i));
end
xx0=zeros(n,1);
lr0=[0 lr];
ww=[0 ww];
for i=1:k1  
    ww1(i)=disturbance_2(ss(i));
end
for i=1:k1
    lr1(i)=ref_signal_2(ss(i));
end
xx00=zeros(n,llr/h);
lr1=[zeros(1,llr/h) lr1];
ww1=[zeros(1,llr/h) ww1];
[yy,uu,um,ee]=simulation(k,lr1,ww1,K_e,K_x,xx00);
[yy1,uu1,um1,ee1]=simulation1(k,lr0,ww,F_e,F_x,F_ksai,Ec,xx0);
ShuChu(k,yy,uu,um,ee,yy1,uu1,um1,ee1,lr0,ww,1);
%% 预见步长为0.60
function [yy,uu,um,ee]=simulation(k,lr,ww,K_e,K_x,xx0)
global A;global B;global C;global Bd;global n; global llr1
global p;global llr;global h;global s;global P; global u0;
xx=zeros(n,k);
xx=[xx0 xx];
uu=zeros(1,k+llr/h);
um=zeros(1,k+llr/h);
yy=C*xx;
lr(:,1:k+llr/h);
ee=zeros(1,k+llr/h);
int_e=zeros(p,k+llr/h);
int_r=zeros(p,k+llr/h);
qq=zeros(1,k+llr/h);
for i=llr/h:k+llr/h-1
    if uu(:,i)>=u0
        um(:,i)=u0;
    elseif uu(:,i)<=-u0
        um(:,i)=-u0;
    else
        um(:,i)=uu(:,i);
    end
    xx(:,i+1)=(h*A+eye(n))*xx(:,i)+h*B*um(:,i)+h*Bd*ww(:,i);
    yy(:,i+1)=C*xx(:,i+1);
    ee(:,i+1)=yy(:,i+1)-lr(:,i+1);
    int_e(:,i+1)=int_e(:,i)+ee(:,i+1);
    int_r(:,i)=0;
    for j=i:i+llr/h
        int_r(:,j+1)=int_r(:,j)+lr(:,j+1);
    end 
    qq(:,i+1)=h*(int_e(:,i+1)-int_r(:,i+llr/h));
    uu(:,i+1)=K_x*xx(:,i+1)-h*K_e*int_r(:,i+llr/h)+h*K_e*int_e(:,i+1);
end

%% 对比
function [yy1,uu1,um1,ee1]=simulation1(k,lr,ww,F_e,F_x,F_ksai,Ec,xx0)
global A;global Bd;global B;global C;global n;global p; global u0;
global h;global M1; global M2;
xx1=zeros(n,k);
xx1=[xx0 xx1];
uu1=zeros(1,k);
um1=zeros(1,k);
yy1=C*xx1;
lr(1,1:k);
ee1=zeros(1,k);
int_e1=zeros(p,k);
int_r1=zeros(p,k);
ksai=zeros(p,k);
for i=1:k
    if uu1(:,i)>=u0
        um1(:,i)=u0;
    elseif uu1(:,i)<=-u0
        um1(:,i)=-u0;
    else
        um1(:,i)=uu1(:,i);
    end
    ksai(:,i+1)=h*ee1(:,i)+h*Ec*(um1(:,i)-uu1(:,i))+ksai(:,i);
    xx1(:,i+1)=(h*A+eye(n))*xx1(:,i)+h*B*um1(:,i)+h*Bd*ww(:,i);
    yy1(:,i+1)=C*xx1(:,i+1);
    ee1(:,i+1)=yy1(:,i+1)-lr(:,i+1);
    ee1(:,2);
%     int_e1(:,i+1)=int_e1(:,i)+ee1(:,i+1);
%     int_r1(:,i+1)=0;
    uu1(:,i+1)=F_e*ee1(:,i+1)+F_x*M1*xx1(:,i+1)+F_ksai*ksai(:,i+1);
end 

%%
%ShuChu(k,yy,yy1,uu,um,uu1,um1,uu2,um2,ee,yy2,lr0,ww,ee1,ee2,1);
function [U,P,L,H,K_e,K_x]=LMI(k)
global A;global B;global C; global Bd;global n;global m;global p;global q;
%%%%%%构造误差系统
A0=[zeros(p,p) C;zeros(n,p) A];
B0=[zeros(p,m);B];
C0=[eye(p,p) zeros(p,n)];
D0=[-eye(p) zeros(p,q);zeros(n,p)  Bd];
%%%%%%求解线性矩阵不等式
gamma=1.5;q_e=1.5; q_x=0.1; r=0.1;
M=[sqrt(q_e*eye(p)) zeros(p,n);zeros(n,p) sqrt(q_x*eye(n));zeros(m,p) zeros(m,n)];
N=[zeros(p,m);zeros(n,m);sqrt(r*eye(m))];
%%%初始化LMI
setlmis([]);
%%%定义决策变量
U=lmivar(1,[m 1]);
P=lmivar(1,[p+n 1]);
H=lmivar(2,[m n+p]);
L=lmivar(2,[m n+p]);
%%%LMI的描述
% 1st LMI
lmiterm([1 1 1 P],A0,1,'s');
lmiterm([1 1 1 L],B0,1,'s');
lmiterm([1 1 3 0],D0);
lmiterm([1 2 1 U],-1,B0');
lmiterm([1 2 1 H],1,1);
lmiterm([1 2 1 L],1,1);
lmiterm([1 2 2 U],-2,1);
lmiterm ([1 3 3 0],-gamma*gamma);
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
% 3rd-4th LMI
lmiterm([-4 1 1 P],1,1);
lmiterm([-5 1 1 U],1,1);
%%%LMI的求解
lmis=getlmis;
[tmin,xfeas]=feasp(lmis);
U=dec2mat(lmis,xfeas,U);
P=dec2mat(lmis,xfeas,P);
H=dec2mat(lmis,xfeas,H);
L=dec2mat(lmis,xfeas,L);
%控制增益阵的求解;
K=L*inv(P);
[k11,k12]=size(K);
K_e=K(:,1:q);
K_x=K(:,q+1:k12);

function [F_x,F_e,F_ksai,Ec]=LMI1(k)
global A;global B;global C; global D;global n;global m;global p;
global q;global u0;global M2;
Ec=zeros(1,m);
%%%%%%构造误差系统
E=[eye(p);zeros(n-p,p)];
A0=[M2*A*inv(M2) zeros(n,p);E' zeros(p,p)];
B01=[M2*B;zeros(p,m)];
B02=[zeros(n,m);Ec];
V=[zeros(n,p);eye(p)];
B03=[eye(n);zeros(p,n)];
%D0=[-eye(p) zeros(p,q);zeros(n,p)  D];
%%%%%%求解线性矩阵不等式
lamda=1.25;
% gamma=1.5;q_e=1.5; q_x=0.1; r=0.1;
% M=[sqrt(q_e*eye(p)) zeros(p,n);zeros(n,p) sqrt(q_x*eye(n));zeros(m,p) zeros(m,n)];
% N=[zeros(p,m);zeros(n,m);sqrt(r*eye(m))];
%%%初始化LMI
setlmis([]);
%%%定义决策变量
Y=lmivar(2,[m n+p]);
W=lmivar(1,[p+n 1]);
M=lmivar(2,[p m]);
L=lmivar(1,[m 1]);
R=lmivar(1,[n 1]);
X=lmivar(2,[n+p m]);
%%%LMI的描述
% 1st LMI
%Lamda=W*A0'+A0*W+B01*Y+Y'*B01'+lamda*W;
lmiterm([1 1 1 W],A0,1,'s');
lmiterm([1 1 1 W],lamda,1);
lmiterm([1 1 1 Y],B01,1,'s');
lmiterm([1 1 2 L],B01,-1);
lmiterm([1 1 2 M],V,-1);
lmiterm([1 1 2 X],1,1);
lmiterm([1 1 3 0],B03);
lmiterm([1 2 2 L],-2,1);
lmiterm ([1 3 3 R],-lamda,1);
%2rd LMI
for i=1:m
    lmiterm([-2 1 1 W],1,1);
    lmiterm([-2 2 1 Y(i,:)],1,1);
    lmiterm([-2 2 1 -X(:,i)],-1,1);
    lmiterm([-2 2 2 0],u0(i,:)*u0(i,:));
end
% 3rd LMI
for i=1:m
    lmiterm([-3 1 1 W],-B03*B03'*A0,1,'s');
    lmiterm([-3 1 1 Y],B03*B03'*B01,1,'s');
    lmiterm([-3 1 2 0],B03);
    lmiterm([-3 3 1 Y(i,:)],1,1);
    lmiterm([-3 2 2 R],lamda,1);
    lmiterm ([-3 3 3 0],(1/lamda)*u0(i,:)*u0(i,:));
end
%4th-6th
lmiterm([-4 1 1 W],1,1);
lmiterm([-5 1 1 R],1,1);
lmiterm([-6 1 1 L],1,1);
%%%LMI的求解
lmis=getlmis;
[tmin,xfeas]=feasp(lmis);
W=dec2mat(lmis,xfeas,W);
Y=dec2mat(lmis,xfeas,Y);
M=dec2mat(lmis,xfeas,M);
L=dec2mat(lmis,xfeas,L);
%控制增益阵的求解;
F=Y*inv(W);
Ec=M*inv(L)
[k11,k12]=size(F);
F_e=F(:,1:q)
F_x=F(:,q+1:q+(n-p))
F_ksai=F(:,q+(n-p)+1:k12)

function ShuChu(k,yy,uu,um,ee,yy1,uu1,um1,ee1,lr,ww,FLG)
global s;global llr1;global h;global ss;global k1;global tf;global llr;
temp=[yy1,ww];MMx=max(temp);MI=min(temp);
figure(FLG);
s=0:h:tf;
    plot(s,lr(:,1:k),'k','LineWidth',1);
    hold on
    plot(s,ww(:,1:k),'k--','LineWidth',1);
    hold on
    plot(s,yy1(:,1:k),'b:','LineWidth',1);
%     hold on
%     plot(s,yy2(:,llr/h+1:k+llr/h-1),'b-.','LineWidth',1);
%     hold on
     plot(s,yy(:,llr/h+1:k+llr/h),'r--','LineWidth',1);
    hold off
   %title('Previewable length：0.60');%修改显示参数
    legend('r(t)','w(t)','y(t): Method in 21','y(t): Preview length: 1.00s')
    xlabel('t/s');
    ylabel('r(t), w(t), y(t)');  
hold off;
ylim([-0.5 MMx+0.1*abs(MMx)])

figure(FLG+1);
%s=0:h:30-h;
    subplot(211)
    temp3=[um1,uu1];MMx3=max(temp3);MI3=min(temp3);
    ylim([MI3-0.02 MMx3+0.1*abs(MMx3)])
    plot(s,uu1(:,1:k),'b-','LineWidth',1);
    hold on
    plot(s,um1(:,1:k),'r--','LineWidth',1);
    hold off
   title('Method in 21');%修改显示参数
    xlabel('t/s');
    ylabel(' u(t) , sat(u(t))'); 
%     subplot(312)
%     temp2=[um2,uu2];MMx2=max(temp2);MI2=min(temp2);
%     plot(s,uu2(:,llr/h+1:k+llr/h-1),'b-','LineWidth',1);
%     hold on
%     plot(s,um2(:,llr/h+1:k+llr/h-1),'r--','LineWidth',1);
%     hold off
%    title('Previewable length：0.48');%修改显示参数    subplot(311)
%     xlabel('t/s');
%     ylabel(' u(t) , sat(u(t))'); 
   subplot(212)
   temp1=[um,uu,uu1,um1];MMx1=max(temp1);MI1=min(temp1);
    plot(s,uu(:,llr/h+1:k+llr/h),'b-','LineWidth',1);
    hold on
    plot(s,um(:,llr/h+1:k+llr/h),'r--','LineWidth',1);
    hold off
   title('Previewable length：1.00s');%修改显示参数
    legend('u(t) ', 'sat(u(t))');
    xlabel('t/s');
    ylabel(' u(t) , sat(u(t))');  
% hold off;
%ylim([MI1-0.1*abs(MI1) MMx1+0.1*abs(MMx1)])
% 
%误差信号图像
MMx=max([ee1,ee]);MI=min([ee1,ee]);
figure(FLG+2)
    plot(s,ee1(:,1:k),'b:','LineWidth',1);
%     hold on
%     plot(s,ee2(:,1+llr1/h:k+llr1/h-1),'b-.','LineWidth',1);
    hold on
    plot(s,ee(:,1+llr/h:k+llr/h),'r--','LineWidth',1);
    hold off
    %title('Mean value of tracking error:',er);
    xlabel('t/s')
    ylabel('tracking error: e(t)')
    hold off;
    legend('Method in 21','Preview length:1.00s')

ylim([MI-0.1*abs(MI) MMx+0.1*abs(MMx)])

 %目标信号1
    function lr1=ref_signal_1(t)     
    lr1=0.5*(t>=0&t<15)+(-0.5).*(t>=15);
    %目标信号2
function lr2=ref_signal_2(t)
    lr2=0.*(t>=0&t<5)+(t-5).*(t>=5&t<10)+5.*(t>=10); 
%干扰信号1  
function ww1=disturbance_1(t)
  ww1=0.1;  
  %干扰信号2
function ww2=disturbance_2(t)
  ww2=0.*(t>=0&t<25)+3.5*(t>=25&t<30);  

