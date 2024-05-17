%1.最优静态工作点理论(画图全英文标注)

clear;
clc
close all;
set(0,'defaultfigurecolor','w')%显示背景设置为白色


Isc=8;% 短路电流,单位：毫安
Vcc=12;%直流电压,单位：伏特
Uoc=Vcc;%开路电压,等于Vcc,单位：伏特
beta=200;%电流放大倍数
Ube0=0.7;%基射正偏电压，单位：伏特
Rb=300000;%300k欧姆；
Rc=Uoc/Isc*1000;%4k欧姆；
Ib0=(Uoc-Ube0)/Rb;%基极电流,单位：安
Ic0=beta*Ib0;%集电极电流,单位：安
Uce0=Vcc-Ic0*Rc;%集射级电压,单位：安
Q0=[Ib0,Ic0,Uce0]%参考Q点,由具体电路参数决定

Uce=0:1/2:Uoc;%基射电压取值范围
Ic=-Isc/Uoc*Uce+Isc;%直流负载线
Ib=Ic/beta;
Ube=Vcc-Ic/beta*Rb;
P=Uce.*Ic;%功率
m=size(Uce,2);
Rc1=Rc-1000*floor(m/2):1000:Rc+1000*floor(m/2);
% Rc2=Rc-1000*floor(m/2):300:Rc+1000*floor(m/2);
P1=Uce.*Uce/4./Rc1;
% P2=Vcc.*Vcc/4./Rc2;

Uth=zeros(1,m);%电压裕度
Un=zeros(1,m);%噪声裕度Vcc/2-Uth
S=zeros(1,m);%电路工作状态
Udev=zeros(1,m);%偏离度

Uce_th=zeros(1,m);%允许Uce相对变化量；
Vcc_th=zeros(1,m);%允许Vcc相对变化量；
Ib_th=zeros(1,m);%允许Vcc相对变化量
Ic_th=zeros(1,m);%允许Vcc相对变化量
Ube_th=zeros(1,m);%允许Vbe相对变化量；
Rb_th=zeros(1,m);%允许Rb相对变化量；
Rc_th=zeros(1,m);%允许Rc相对变化量；
beta_th=zeros(1,m);%允许beta相对变化量；
small=1.0e-6*rand(1);%加一个小量避免分母为零，产生奇异
for i=1:m
    S(:,i)=Uce(:,i)/(Uoc/2);
    if Uce(:,i)<=Uoc/2
        Uth(:,i)=Uce(:,i);
        Udev(:,i)=1-S(:,i);
    else
        Uth(:,i)=Uoc-Uce(:,i);
        Udev(:,i)=S(:,i)-1;
    end
    Un(:,i)=Vcc/2-Uth(:,i);
    Uce_th(:,i)=1-Udev(:,i);
    temp=Vcc-Uce(:,i)*(1+Uce_th(:,i));
    Ube_th(:,i)=(Vcc-abs(temp)*Rb/Rc/beta)/(Ube(:,i))-1;
    Vcc_th(:,i)=(Ic(:,i)*Rc/1000+Uce(:,i)*(1+Uce_th(:,i)))/Vcc-1;
    Ib_th(:,i)=abs(temp)/Rc/beta/Ib(:,i)*1000-1;%乘以1000是量纲对齐，注意量纲
    Ic_th(:,i)=abs(temp)/Rc/Ic(:,i)*1000-1;%乘以1000是量纲对齐，注意量纲
    Rb_th(:,i)=Rc*beta*(Vcc-Ube(:,i))/abs(temp)/Rb-1;
    Rc_th(:,i)=abs(temp)/Ic(:,i)*1000/Rc-1;
    beta_th(:,i)=abs(temp)/Rc/beta/Ib(:,i)*1000-1;%乘以1000是量纲对齐，注意量
end


figure(1)
subplot(2,1,1);
colororder({'r','b'})
yyaxis left
plot(Uce,Ic,'r^-','LineWidth',1.2);
xlabel('Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('Ie(mA)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');

yyaxis right
hold on;
plot(Uce,Uth,'b>-','LineWidth',1.2);
plot(Uce,Un,'b*-','LineWidth',1.2);
ylabel('Uth(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
str={'DC Load Line:Ic=-Isc/Uoc*Uce+Isc','Voltage Margin Line:Uth','Noise Margin Line:Vcc/2-Uth'};
legend(str);
grid on;

subplot(2,1,2);
colororder({'r','b'})
yyaxis left
plot(Uce,Ic,'r^-','LineWidth',1.2);
xlabel('Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('Ie(mA)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');

yyaxis right
hold on
plot(Uce,S,'b<-','LineWidth',1.2);
plot(Uce,Udev,'bo-','LineWidth',1.2);
str={'DC Load Line:Ic=-Isc/Uoc*Uce+Isc','Circuit Operating Status:S','Skewness:Udev'};
legend(str);
grid on;


figure(2)
colororder({'r','b'})
yyaxis left
plot(Uce,Ic,'r^-','LineWidth',1.2);
xlabel('Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('Ie(mA)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');

yyaxis right
hold on;
plot(Uce,P,'b*-','LineWidth',1.2);
ylabel('P(W)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
str={'DC Load Line:Ic=-Isc/Uoc*Uce+Isc','Power:P=Uce*Ic'};
legend(str);
grid on;

%创建第三套坐标轴，设定颜色及位置
ax2=axes('Position',get(gca,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','m','YColor','m');
hold on;
plot(Rc1,P1,'m*-','LineWidth',1.2);
xlabel('Rc(\Omega)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
str={'Power:P=Uce*Uce/(4*Rc)'};
legend(str);
grid on;

figure(3)
subplot(4,2,1)
plot(Uce,Uce_th,'r.-','LineWidth',1.2);
xlabel('(a) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\deltaUce','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;

subplot(4,2,2)
plot(Uce,Ube_th,'g>-','LineWidth',1.2);
xlabel('(b) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\deltaUbe','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;

subplot(4,2,3)
plot(Uce,Vcc_th,'b^-','LineWidth',1.2);
xlabel('(c) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\deltaVcc','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;

subplot(4,2,4)
plot(Uce,Rb_th,'c*-','LineWidth',1.2);
xlabel('(d) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\deltaRb','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;

subplot(4,2,5)
plot(Uce,Ic_th,'r<-','LineWidth',1.2);
xlabel('(e) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\deltaIc','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;

subplot(4,2,6)
plot(Uce,Ib_th,'g^-','LineWidth',1.2);
xlabel('(f) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\deltaIb','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;

subplot(4,2,7)
plot(Uce,Rc_th,'b^-','LineWidth',1.2);
xlabel('(g) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\deltaRc','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;

subplot(4,2,8)
plot(Uce,beta_th,'c^-','LineWidth',1.2);
xlabel('(h) Uce(V)','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
ylabel('\delta\beta','FontSize',12,'Fontname', 'times new roman','FontWeight','bold');
grid on;
