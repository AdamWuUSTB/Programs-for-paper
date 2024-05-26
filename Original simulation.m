function SPC
clear;clc;
global A;global D;global B;global C;global tf;global s;global D1;global B1;global C1;
global n;global m;global p;global q;global h;global llr;global llr1;global u0;
llr=1.00;llr1=0.50;%设定预见长度
% A=[-0.8 -0.0006 -12 0;0 -0.014 -16.64 -32.2;1 0.0001 -1.5 0;1 0 0 0];
% B=[-19 -3;-0.66 -0.5;-0.16 -0.5;0 0];C=[0 0 0 1;0 0 -1 1]; D=[1 -3;0 -5;0.16 0.5;1 2];
A=[-2.98 0.4 0 -0.034;-0.99 -0.21 0.035 -0.0011;0 0 0 1;0.39 -5.55 0 -1.89];
B=[-0.032;0;0;-1.6];D=[0.5;0;0.6;0];%飞行器参数取0.4
C=[0 0 1 0];
u0=5;
[n,n]=size(A); [n,m]=size(B); [p,n]=size(C); [n,q]=size(D);
t0=0; tf=50; h=0.001; s=t0:h:tf; [k0,k]=size(s);ss=t0:h:tf+llr+llr;[k0,k1]=size(ss);
[U,P,L,H,K_e,K_x] = LMI(k);
for i=1:k1  
    ww(i)=disturbance_2(ss(i));
end
for i=1:k1
    lr(i)=ref_signal_2(ss(i));
end
xx0=zeros(n,llr/h);
lr0=[zeros(1,llr/h) lr];
ww=[zeros(1,llr/h) ww];
[yy,uu,um,yy1,yy2,uu1,uu2,um1,um2,ee,ee1,ee2]=simulation(k,lr0,ww,K_e,K_x,xx0);
ShuChu(k,yy,yy1,uu,um,uu1,um1,uu2,um2,ee,yy2,lr0,ww,ee1,ee2,1);
%% 预见步长为0.60
function [yy,uu,um,yy1,yy2,uu1,uu2,um1,um2,ee,ee1,ee2]=simulation(k,lr,ww,K_e,K_x,xx0)
global A;global B;global C;global D;global n; global llr1
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
    xx(:,i+1)=(h*A+eye(n))*xx(:,i)+h*B*um(:,i)+h*D*ww(:,i);
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
%% 不含可预见参考信息-预见步长为0
xx1=zeros(n,k);
xx1=[xx0 xx1];
uu1=zeros(1,k+llr/h);
um1=zeros(1,k+llr/h);
yy1=C*xx1;
lr(:,1:k+llr/h);
ee1=zeros(1,k+llr/h);
int_e1=zeros(p,k+llr/h);
int_r1=zeros(p,k+llr/h);
for i=llr/h:k+llr/h-1
    if uu1(:,i)>=u0
        um1(:,i)=u0;
    elseif uu1(:,i)<=-u0
        um1(:,i)=-u0;
    else
        um1(:,i)=uu1(:,i);
    end
    xx1(:,i+1)=(h*A+eye(n))*xx1(:,i)+h*B*um1(:,i)+h*D*ww(:,i);
    yy1(:,i+1)=C*xx1(:,i+1);
    ee1(:,i+1)=yy1(:,i+1)-lr(:,i+1);
    int_e1(:,i+1)=int_e1(:,i)+ee1(:,i+1);
    int_r1(:,i+1)=0;
    uu1(:,i+1)=K_x*xx1(:,i+1)+h*K_e*int_e1(:,i+1)-h*K_e*int_r(:,i+1);
end
%% 预见步长为0.42
xx2=zeros(n,k);
xx2=[xx0 xx2];
uu2=zeros(1,k+llr1/h);
um2=zeros(1,k+llr1/h);
yy2=C*xx2;
lr(:,1:k+llr1/h);
ee2=zeros(1,k+llr1/h);
int_e2=zeros(p,k+llr1/h);
int_r2=zeros(p,k+llr1/h);
qq=zeros(1,k+llr1/h);
for i=llr/h:k+llr/h-1
    if uu2(:,i)>=u0
        um2(:,i)=u0;
    elseif uu2(:,i)<=-u0
        um2(:,i)=-u0;
    else
        um2(:,i)=uu2(:,i);
    end
    xx2(:,i+1)=(h*A+eye(n))*xx2(:,i)+h*B*um2(:,i)+h*D*ww(:,i);
    yy2(:,i+1)=C*xx2(:,i+1);
    ee2(:,i+1)=yy2(:,i+1)-lr(:,i+1);
    int_e2(:,i+1)=int_e2(:,i)+ee2(:,i+1);
    int_r2(:,i)=0;
    for j=i:i+llr1/h
        int_r2(:,j+1)=int_r2(:,j)+lr(:,j+1);
    end 
    qq2(:,i+1)=h*(int_e2(:,i+1)-int_r2(:,i+llr1/h));
    uu2(:,i+1)=K_x*xx2(:,i+1)-h*K_e*int_r2(:,i+llr1/h)+h*K_e*int_e2(:,i+1);
end
%%
%ShuChu(k,yy,yy1,uu,um,uu1,um1,uu2,um2,ee,yy2,lr0,ww,ee1,ee2,1);

function [U,P,L,H,K_e,K_x]=LMI(k)
global A;global B;global C; global D;global n;global m;global p;global q;
%%%%%%构造误差系统
A0=[zeros(p,p) C;zeros(n,p) A];
B0=[zeros(p,m);B];
C0=[eye(p,p) zeros(p,n)];
D0=[-eye(p) zeros(p,q);zeros(n,p)  D];
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
P
L
K=L*inv(P)
[k11,k12]=size(K);
K_e=K(:,1:q)
K_x=K(:,q+1:k12)

function ShuChu(k,yy,yy1,uu,um,uu1,um1,uu2,um2,ee,yy2,lr,ww,ee1,ee2,FLG)
global s;global llr1;global h;global ss;global k1;global tf;global llr;
temp=[yy1,yy,ww];MMx=max(temp);MI=min(temp);
figure(FLG);
s=0:h:tf-h;
    plot(s,lr(:,llr/h+1:k+llr/h-1),'k','LineWidth',1);
    hold on
    plot(s,ww(:,llr/h+1:k+llr/h-1),'k--','LineWidth',1);
    hold on
    plot(s,yy1(:,llr/h+1:k+llr/h-1),'r:','LineWidth',1);
    hold on
    plot(s,yy2(:,llr/h+1:k+llr/h-1),'b-.','LineWidth',1);
    hold on
    plot(s,yy(:,llr/h+1:k+llr/h-1),'g--','LineWidth',1);
    hold off
   %title('Previewable length：0.60');%修改显示参数
    legend('r(t)','w(t)','y(t):0.00s','y(t):0.50s','y(t):1.00s')
    xlabel('t/s');
    ylabel('r(t), w(t), y(t)');  
hold off;
ylim([-0.5 MMx+0.1*abs(MMx)])

figure(FLG+1);

%s=0:h:30-h;
    subplot(311)
    temp3=[um2,uu2];MMx3=max(temp3);MI3=min(temp3);
    ylim([MI3-2 MMx3+0.1*abs(MMx3)])
    plot(s,uu1(:,llr/h+1:k+llr/h-1),'b-','LineWidth',1);
    hold on
    plot(s,um1(:,llr/h+1:k+llr/h-1),'r--','LineWidth',1);
    hold off
   title('Previewable length：0.00s');%修改显示参数
    xlabel('t/s');
    ylabel(' u(t) , sat(u(t))'); 
    subplot(312)
    temp2=[um2,uu2];MMx2=max(temp2);MI2=min(temp2);
    plot(s,uu2(:,llr/h+1:k+llr/h-1),'b-','LineWidth',1);
    hold on
    plot(s,um2(:,llr/h+1:k+llr/h-1),'r--','LineWidth',1);
    hold off
   title('Previewable length：0.50s');%修改显示参数    subplot(311)
    xlabel('t/s');
    ylabel(' u(t) , sat(u(t))'); 
   subplot(313)
   temp1=[um,uu];MMx1=max(temp1);MI1=min(temp1);
    plot(s,uu(:,llr/h+1:k+llr/h-1),'b-','LineWidth',1);
    hold on
    plot(s,um(:,llr/h+1:k+llr/h-1),'r--','LineWidth',1);
    hold off
   title('Previewable length：1.00s');%修改显示参数
   legend('u(t) ', 'sat(u(t))');
    xlabel('t/s');
    ylabel(' u(t) , sat(u(t))');  
hold off;
ylim([MI1-0.1*abs(MI1) MMx1+0.1*abs(MMx1)])
% 
%误差信号图像
MMx=max([ee,ee1,ee2]);MI=min([ee,ee1,ee2]);
figure(FLG+2)
    plot(s,ee1(:,1+llr1/h:k+llr1/h-1),'r:','LineWidth',1);
    hold on
    plot(s,ee2(:,1+llr1/h:k+llr1/h-1),'b-.','LineWidth',1);
    hold on
    plot(s,ee(:,1+llr1/h:k+llr1/h-1),'g--','LineWidth',1);
    hold off
    %title('Mean value of tracking error:',er);
    xlabel('t/s')
    ylabel('tracking error: e(t)')
    hold off;
    legend('e(t):0.00s','e(t):0.50s','e(t):1.00s')

ylim([MI-0.1*abs(MI) MMx+0.1*abs(MMx)])

 %目标信号1
    function lr1=ref_signal_1(t)     
    lr1=1.*sin(t)+2;
    %目标信号2
function lr2=ref_signal_2(t)
    lr2=0.*(t>=0&t<5)+(t-5).*(t>=5&t<10)+5.*(t>=10); 
%干扰信号1  
function ww1=disturbance_1(t)
  ww1=0.*(t>=0&t<7)+3.5*(t>=7&t<10);   
  %干扰信号2
function ww2=disturbance_2(t)
  ww2=0.*(t>=0&t<25)+3.5*(t>=25&t<30);   

