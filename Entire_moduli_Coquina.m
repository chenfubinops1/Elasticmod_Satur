clc
clear
load('Ksatp.mat');
load('Gsatp.mat');
K_d = Ksatp;
mu = Gsatp;
[nf,ns] = size(K_d);
%% plot
MkSize=12;
LdWidth=2;
FtSize=20;
FtName='Arial';
FtSize2=16;
fm=[200 400 700 1000 2000 4000 7000 10000 20000 40000 70000 100000];
Kgly_P1=[28.219 29.003 29.357 29.445 30.267 31.386 31.99 32.898 35.212 37.972 40.405 42.402]*1e9;
Kgly_P2=[32.088 33.624 34.961 35.419 35.786 36.647 37.862 39.248 40.246 43.297 46.357 47.983]*1e9;
Kgly_P4=[37.462 39.262 41.448 42.631 43.728 45.728 46.971 48.299 52.141 55.121 56.997 57.271]*1e9;
QK1=[0.139 0.145 0.138 0.127 0.117 0.103 0.159 0.15 0.177 0.141 0.078 0.018];
QK2=[0.139 0.145 0.138 0.127 0.117 0.103 0.159 0.15 0.177 0.141 0.078 0.018];
QK3=[0.055 0.041 0.053 0.045 0.085 0.115 0.172 0.195 0.182 0.124 0.075 0.063];

Kdry=[7.38 10.25 12.7 15.65 18.32]'*1e9;


yita = 1;  %%Glycerine
k0 = 5*10^-17;
% K_d = 4*10^9;
K_s = 77e9;
K_f = 4.36e9;
% mu = 19.6*10^9;
por = 0.075;
rho_s = 2540;
rho_f = 1261;
rho12 = -83;
rho11 = (1-por)*rho_s-rho12;
rho22 = por*rho_f-rho12;
rho_total = rho_s*(1-por)+rho_f*por;


b = yita*por^2/k0; %%viscosity


%% Gassmann

K_gassm26 = K_d(1,26)+(1-K_d(1,26)/K_s)^2./(por/K_f+(1-por)/K_s-K_d(1,26)/K_s^2);
%
% v0 = sqrt((K_sat+4*mu/3)/(2650*(1-f)+1000*f));

n = 1;
for nn = 1:100:401
    %%% 应该考虑骨架应力对骨架弹性参数的影响
    for m = 1:nf
        
        w = 10^((m-41)/10);
        
        tao = 2;
        k_w = k0/((1-0.5i*tao*k0*rho_f*w/yita/por)^(0.5)-1i*tao*k0*rho_f*w/yita/por); %% 动态渗透率
        theta = 1i*k_w/yita/w;
        arf = 1-K_d(m,nn)/K_s;
        bita = (arf-por)/K_s+por/K_f;
        bs = rho_f*theta*w^2;
        b0 = -bita*(K_d(m,nn)+4*mu(m,nn)/3+arf^2/bita)/arf;
        c = (arf-bs*rho_total/(rho_f*b0))/(arf+bs);
        K_gassm =  K_d(m,nn)+(1-K_d(m,nn)/K_s)^2./(por/K_f+(1-por)/K_s-K_d(m,nn)/K_s^2);
        %         K_gassm = K_d(m,nn)+arf^2/((arf-por)/K_s+por/K_f);
        b_plus = 0.5*b0*(c-sqrt(c^2-4*arf*(1-c)/b0));
        b_sub = 0.5*b0*(c+sqrt(c^2-4*arf*(1-c)/b0));
        rho_w = rho_total+rho_f^2*w^2*theta;  %% w应该为w_f
        
        kp0 = w/sqrt((K_gassm+4*mu(m,nn)/3)/rho_total);
        
        kp1(n,m) = kp0*sqrt((1+b_plus*rho_f/rho_total)/(1-b_plus/b0));
        kp2(n,m) = kp0*sqrt((1+b_sub*rho_f/rho_total)/(1-b_sub/b0));
        ks(n,m) = w*sqrt(rho_w/mu(m,nn));
        
        vp1(n,m) = w/real(kp1(n,m));
        vp2(n,m) = w/real(kp2(n,m));
        vs(n,m) = w/real(ks(n,m));
        
        Q1(n,m) = imag(kp1(n,m))/(real(kp1(n,m)));
        Q2(n,m) = imag(kp2(n,m))/(real(kp2(n,m)));
        Qs(n,m) = imag(ks(n,m))/(real(ks(n,m)));
        
        %% 体积模量
        vp1_com(n,m) = w/kp1(n,m);
        vs_com(n,m) = w/ks(n,m);
        
        K_satp(n,m) = rho_total*(vp1_com(n,m)^2-(4/3)*vs_com(n,m)^2);
        Mu_satp(n,m) = rho_total*vs(n,m)^2;
        
        QK(n,m) = imag(K_satp(n,m))/(real(K_satp(n,m)));
    end
    
    n = n+1;
end


figure
%%fast Vp
y = -6:0.1:10;

figure(1)
set(gcf,'Position',[500 150 800 600])
plot(y,vp1(1,:)/1e3,'color',[0 0 0],'linewidth',LdWidth)
hold on
plot(y,vp1(2,:)/1e3,'color',[0 0 1],'linewidth',LdWidth)
hold on
plot(y,vp1(3,:)/1e3,'color',[0.196 0.803 0.196],'linewidth',LdWidth)
hold on
plot(y,vp1(5,:)/1e3,'color',[1 0 0],'linewidth',LdWidth)
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('Fast P-wave velocity (km/s)','FontName',FtName,'FontSize',FtSize);
axis([-2 10 4 6.5])
set(gca,'xtick',-4:2:10,'ytick',4:0.5:6.5,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa','Location','Southeast','Numcolumns',2);%'Numcolumns',2
legend('boxoff')
title(lgd,'Coquina limestone \phi = 7.5%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')


%%slow Vp
figure(2)
set(gcf,'Position',[500 150 800 600])
plot(y,vp2(1,:)/1e3,'color',[0 0 0],'linewidth',LdWidth)
hold on
plot(y,vp2(2,:)/1e3,'color',[0 0 1],'linewidth',LdWidth)
hold on
plot(y,vp2(3,:)/1e3,'color',[0.196 0.803 0.196],'linewidth',LdWidth)
hold on
plot(y,vp2(5,:)/1e3,'color',[1 0 0],'linewidth',LdWidth)
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('Slow P-wave velocity (km/s)','FontName',FtName,'FontSize',FtSize);
axis([-2 10 0 1.25])
set(gca,'xtick',-4:2:10,'ytick',0:0.25:1.25,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa','Location','Southeast');%'Numcolumns',2
legend('boxoff')
title(lgd,'Coquina limestone \phi = 7.5%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')

%% Vs
figure(3)
semilogx(y,vs(1,:),'color',[0 0 0],'linewidth',1.7);
hold on
semilogx(y,vs(2,:),'color',[0 0 1],'linewidth',1.7);
hold on
semilogx(y,vs(3,:),'color',[0.196 0.803 0.196],'linewidth',1.7);
hold on
semilogx(y,vs(4,:),'color',[1 0 0],'linewidth',1.7);
legend('1','2','3','4');
xlabel('Frequency/Hz');
ylabel('Vs/m·s(-1)');
axis([0 10^9 -inf inf])
%% K_sat
figure(4)
set(gcf,'Position',[500 150 800 600])
plot(y,K_satp(1,:)/(1e9),'color',[0 0 0],'linewidth',LdWidth)
hold on
plot(y,K_satp(2,:)/(1e9),'color',[0 0 1],'linewidth',LdWidth)
hold on
plot(y,K_satp(3,:)/(1e9),'color',[1 0 0],'linewidth',LdWidth)
hold on
plot(log10(fm/1e3),Kgly_P4/1e9,'o','MarkerSize',MkSize,'MarkerEdgecolor',[0 0 1],'MarkerFacecolor',[0 0 1]);
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('Bulk modulus (GPa)','FontName',FtName,'FontSize',FtSize);
axis([-2 4 20 60])
set(gca,'xtick',-2:1:4,'ytick',20:10:60,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,  10 MPa',' Model,  20 MPa',' Experimental,  10 MPa','Location','Southeast');%'Numcolumns',2
legend('boxoff')
title(lgd,'Coquina limestone \phi = 7.5%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')

%% Mu
figure(5)
semilogx(y,Mu_satp(1,:)/(1e9),'color',[0 0 0],'linewidth',1.7);
hold on
semilogx(y,Mu_satp(2,:)/(1e9),'color',[0 0 1],'linewidth',1.7);
hold on
semilogx(y,Mu_satp(3,:)/(1e9),'color',[0.196 0.803 0.196],'linewidth',1.7);
hold on
semilogx(y,Mu_satp(4,:)/(1e9),'color',[1 0 0],'linewidth',1.7);
legend('1','2','3','4');
xlabel('Frequency/Hz');
ylabel('Mu (MPa)');
axis([0 10^9 -inf inf])



%%%%%%%% ################################
%% fast Vp1 1/Q
figure(6)
set(gcf,'Position',[500 150 800 600])
plot(y,abs(Q1(1,:)),'color',[0 0 0],'linewidth',LdWidth)
hold on
plot(y,Q1(2,:),'color',[0 0 1],'linewidth',LdWidth)
hold on
plot(y,Q1(3,:),'color',[0.196 0.803 0.196],'linewidth',LdWidth)
hold on
plot(y,Q1(4,:),'color',[1 0 0],'linewidth',LdWidth)
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('1/Q_{P1}','FontName',FtName,'FontSize',FtSize);
axis([-2 10 0 0.06])
set(gca,'xtick',-2:2:10,'ytick',0:0.01:0.06,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa','Location','Southeast');%'Numcolumns',2
legend('boxoff')
title(lgd,'Indiana limestone \phi = 13.9%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')

%%slow Vp2 1/Q
figure(7)
set(gcf,'Position',[500 150 800 600])
plot(y,abs(Q2(1,:)),'color',[0 0 0],'linewidth',LdWidth)
hold on
plot(y,abs(Q2(2,:)),'color',[0 0 1],'linewidth',LdWidth)
hold on
plot(y,abs(Q2(3,:)),'color',[0.196 0.803 0.196],'linewidth',LdWidth)
hold on
plot(y,abs(Q2(4,:)),'color',[1 0 0],'linewidth',LdWidth)
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('1/Q_{P2}','FontName',FtName,'FontSize',FtSize);
axis([-2 10 0 1.33])
set(gca,'xtick',-2:2:10,'ytick',0:0.25:1.25,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa','Location','Southeast');%'Numcolumns',2
legend('boxoff')
title(lgd,'Indiana limestone \phi = 13.9%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')



%% Vs 1/Q
figure(8)
semilogx(y,abs(Qs(1,:)),'color',[0 0 0],'linewidth',1.7);
hold on
semilogx(y,abs(Qs(2,:)),'color',[0 0 1],'linewidth',1.7);
hold on
semilogx(y,abs(Qs(3,:)),'color',[0.196 0.803 0.196],'linewidth',1.7);
hold on
semilogx(y,abs(Qs(4,:)),'color',[1 0 0],'linewidth',1.7);
legend('1','2','3','4');
xlabel('Frequency/Hz');
ylabel('1/Qs');
axis([0 10^9 -inf inf])

%% K 1/Q
figure(9)
semilogx(fm/1e3,QK1,'o','MarkerSize',MkSize,'MarkerFacecolor','k','MarkerEdgecolor','k')
hold on
%     semilogx(fm/1e3,QK2,'d','MarkerSize',MkSize,'MarkerFacecolor','r','MarkerEdgecolor','r')
%     hold on
semilogx(fm/1e3,QK3,'s','MarkerSize',MkSize,'MarkerFacecolor','b','MarkerEdgecolor','b')
hold on
semilogx(y,abs(QK(1,:)),'color',[0 0 0],'linewidth',1.7);
hold on
semilogx(y,abs(QK(2,:)),'color',[0 0 1],'linewidth',1.7);
hold on
semilogx(y,abs(QK(3,:)),'color',[0.196 0.803 0.196],'linewidth',1.7);
hold on
semilogx(y,abs(QK(4,:)),'color',[1 0 0],'linewidth',1.7);
legend('1','2','3','4');
xlabel('Frequency/Hz');
ylabel('1/QK');
axis([0 10^9 -inf inf])


figure(10)
set(gcf,'Position',[500 150 800 600])
plot(y,abs(QK(1,:)),'color',[0 0 0],'linewidth',LdWidth)
hold on
plot(y,abs(QK(2,:)),'color',[0 0 1],'linewidth',LdWidth)
hold on
plot(y,abs(QK(3,:)),'color',[1 0 0],'linewidth',LdWidth)
hold on
plot(log10(fm/1e3),QK3,'o','MarkerSize',MkSize,'MarkerEdgecolor',[0 0 1],'MarkerFacecolor',[0 0 1]);
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('1/Q_{K}','FontName',FtName,'FontSize',FtSize);
axis([-2 4 0 0.4])
set(gca,'xtick',-2:1:4,'ytick',0:0.1:0.4,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,  10 MPa',' Model,  20 MPa',' Experimental,  10 MPa','Location','Southeast');%'Numcolumns',2
legend('boxoff')
title(lgd,'Coquina limestone \phi = 7.5%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')
