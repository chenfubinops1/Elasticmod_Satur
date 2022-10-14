clc
clear
load('Ksatp_india.mat');
load('Gsatp_india.mat');
K_d = Ksatp;
mu = Gsatp;
[nf,ns] = size(K_d);
%% plot
MkSize=12;
LdWidth=2;
FtSize=20;
FtName='Arial';
FtSize2=16;
fm_3_5MPa=[0.00697 0.00999 0.01969 0.0396 0.0679 0.09728 0.111878 0.19966 0.39385 0.7033 0.98807];
fm_5MPa=[0.00697 0.00999 0.01969 0.0396 0.0679 0.09728 0.111878 0.19966 0.39385 0.7033 0.98807];
fm_10MPa=[0.00695 0.00994 0.00196 0.0395 0.0691544 0.1073 0.1995 0.3937];
fm_20MPa=[0.006947 0.00955 0.019619 0.03949 0.069 0.097];
K_3_5MPa=[19.56 19.78 21.43 22.5275 24.5055 26.1538 26.7033 30.5495 34.3956 37.5824 39.3407]*1e9;
K_5MPa=[30.769 31.97 34.0659 36.4835 38.2418 38.9011 39.2308 40.329 40.879 40.3297 40]*1e9;
K_10MPa=[37.4725 42.8571 43.0769 42.6374 42.63 43.5165 42.1978 41.4286]*1e9;
K_20MPa=[41.32 44.62 44.18 43.74 43.956 43.516]*1e9;

QK_3_5MPa=[0.1105 0.1526 0.1539 0.2 0.2368 0.2605 0.2578 0.2605 0.215789 0.14868 0.071];
QK_5MPa=[0.125 0.0868 0.118 0.0973 0.0921 0.0789 0.07368 0.03684 -0.021 -0.0789 -0.07894];


yita = 1;  %%Glycerine
k0 = 3e-13;%% 模拟的时候用3e-13； 验证用真实的3e-17
% K_d = 4*10^9;
K_s = 77e9;
K_f = 4.36e9;
% mu = 19.6*10^9;
por = 0.139;
rho_s = 2345;
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
for nn = 1:50:201
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
axis([-4 7 3.5 6.1])
set(gca,'xtick',-6:1:10,'ytick',3.5:0.5:6,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa','Location','Southeast','Numcolumns',2);%'Numcolumns',2
legend('boxoff')
title(lgd,'Indiana limestone \phi = 13.9%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')

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
axis([-4 7 0 1.25])
set(gca,'xtick',-6:1:10,'ytick',0:0.25:1.5,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa','Location','Southeast');%'Numcolumns',2
legend('boxoff')
title(lgd,'Indiana limestone \phi = 13.9%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')

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
plot(y,K_satp(3,:)/(1e9),'color',[0.196 0.803 0.196],'linewidth',LdWidth)
hold on
plot(y,K_satp(5,:)/(1e9),'color',[1 0 0],'linewidth',LdWidth)
hold on
plot(log10(fm_5MPa),K_5MPa/1e9,'o','MarkerSize',MkSize,'MarkerEdgecolor',[0 0 1],'MarkerFacecolor',[0 0 1]);
hold on
plot(log10(fm_10MPa),K_10MPa/1e9,'o','MarkerSize',MkSize,'MarkerEdgecolor',[0.196 0.803 0.196],'MarkerFacecolor',[0.196 0.803 0.196]);
hold on
plot(log10(fm_20MPa),K_20MPa/1e9,'o','MarkerSize',MkSize,'MarkerEdgecolor',[1 0 0],'MarkerFacecolor',[1 0 0]);
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('Bulk modulus (GPa)','FontName',FtName,'FontSize',FtSize);
axis([-4 2 10 50])
set(gca,'xtick',-4:1:2,'ytick',10:10:50,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa',' Experimental,   5  MPa',' Experimental,  10 MPa',' Experimental,  20 MPa','Location','Southeast','Numcolumns',2);%'Numcolumns',2
legend('boxoff')
title(lgd,'Indiana limestone \phi = 13.9%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')

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
axis([-4 7 0 0.06])
set(gca,'xtick',-4:1:6,'ytick',0:0.01:0.06,'FontName',FtName,'FontSize',FtSize)
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
axis([-4 7 0 1.25])
set(gca,'xtick',-4:1:6,'ytick',0:0.25:1.25,'FontName',FtName,'FontSize',FtSize)
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
set(gcf,'Position',[500 150 800 600])
plot(y,abs(QK(1,:)),'color',[0 0 0],'linewidth',LdWidth)
hold on
plot(y,abs(QK(2,:)),'color',[0 0 1],'linewidth',LdWidth)
hold on
plot(y,abs(QK(3,:)),'color',[0.196 0.803 0.196],'linewidth',LdWidth)
hold on
plot(y,abs(QK(5,:)),'color',[1 0 0],'linewidth',LdWidth)
hold on
plot(log10(fm_3_5MPa),QK_3_5MPa,'o','MarkerSize',MkSize,'MarkerEdgecolor',[0.196 0.803 1],'MarkerFacecolor',[0.196 0.803 1]);
hold on
plot(log10(fm_5MPa),QK_5MPa,'o','MarkerSize',MkSize,'MarkerEdgecolor',[0 0 1],'MarkerFacecolor',[0 0 1]);
xlabel('Log[Frequency] (Hz)','FontName',FtName,'FontSize',FtSize);
ylabel('1/Q_{K}','FontName',FtName,'FontSize',FtSize);
axis([-4 2 -0.12 0.48])
set(gca,'xtick',-4:1:2,'ytick',-0.1:0.1:0.45,'FontName',FtName,'FontSize',FtSize)
lgd=legend(' Model,   0  MPa',' Model,   5  MPa',' Model,  10 MPa',' Model,  20 MPa',' Experimental,  3.5 MPa',' Experimental,    5  MPa','Location','Southeast','Numcolumns',2);%'Numcolumns',2
legend('boxoff')
title(lgd,'Indiana limestone \phi = 13.9%','FontName',FtName,'FontSize',FtSize,'FontWeight','Normal')


