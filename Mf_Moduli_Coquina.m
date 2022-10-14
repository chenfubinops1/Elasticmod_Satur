close all;
clear;
clc;

%% plot
MkSize=12;
LdWidth=2;
FtSize=20;
FtName='Arial';
FtSize2=16;

%%% Basic parameters for modelling
Ks=77e9;
por=0.075;   %
etaw=1e-3;
etagly=1;
Kgly=4.36e9;
p=[3.5 5 10 20 25]'*1e6;

fm=[200 400 700 1000 2000 4000 7000 10000 20000 40000 70000 100000];
Kgly_P1=[28.219 29.003 29.357 29.445 30.267 31.386 31.99 32.898 35.212 37.972 40.405 42.402]*1e9;
Kgly_P2=[32.088 33.624 34.961 35.419 35.786 36.647 37.862 39.248 40.246 43.297 46.357 47.983]*1e9;
Kgly_P4=[37.462 39.262 41.448 42.631 43.728 45.728 46.971 48.299 52.141 55.121 56.997 57.271]*1e9;
QK1=[0.139 0.145 0.138 0.127 0.117 0.103 0.159 0.15 0.177 0.141 0.078 0.018];
QK2=[0.139 0.145 0.138 0.127 0.117 0.103 0.159 0.15 0.177 0.141 0.078 0.018];
QK3=[0.055 0.041 0.053 0.045 0.085 0.115 0.172 0.195 0.182 0.124 0.075 0.063];

Kdry=[7.38 10.25 12.7 15.65 18.32]'*1e9;

Kh=54e9;
Gh=40e9;
vh=(3*Kh-2*Gh)/(6*Kh+2*Gh);
Kex=Kh./Kdry-1;
fex=fit(p,Kex,'exp1');
P1=(0:0.1:180)*1e6;                   %% 用于拟合K_dry的压力范围
Kdry_fit=Kh./(1+fex.a*exp(fex.b*P1));  %% 指数拟合干燥K_dry
cc=16/9*(1-vh^2)/(1-2*vh);  
gamma0=fex.a/cc;   %%指数拟合的K_dry代入（11）in Zimmer得到初始压力时裂缝密度
gamma=gamma0*exp(fex.b*P1); %%裂缝密度随压力
ard_initial_fit=3/(4*pi)*cc*P1/Kh;   %% ard-aspect artio星号in(21) in Zimmer
c_initial_fit=4*pi/3*(-fex.b).*gamma.*P1;
ard_min=4/(3*pi)*(1-vh^2)./(Kh*(1-2*vh))*P1;   %% ar星，与ard_initial_fit一致
cuml_MK_fit=cumtrapz(ard_initial_fit,c_initial_fit);

[mP,nP]=size(P1);
for k=1:(nP-3)
    Khost=Kh;
    ard_initial_fitP=ard_initial_fit((ard_initial_fit-ard_min(k))/ard_min(k)>1e-6); 
    c_initialP=c_initial_fit((ard_initial_fit-ard_min(k))/ard_min(k)>1e-6);  
    [mxP,nxP]=size(c_initialP);        
    ardP=ard_initial_fitP-ard_min(k);  
    
    f=10.^(-6:0.1:10); 
    [mf,nf]=size(f);
    Gfwp=1i*2*pi*f*etagly;
    KfwP = fluid_mod( Kgly, Gfwp', ardP);
    [pf,qf]=Wu(Kh,Gh,KfwP,0,ardP);    
    SpKf=c_initialP.*(1-ard_min(k)./ard_initial_fitP).*(1-KfwP./Kh).*pf;
%     SpKf=SpKf1;
    SpGf=c_initialP.*(1-ard_min(k)./ard_initial_fitP).*qf;
%     SpGf=SpGf1  
    for j=1:nf
    intafK=cumtrapz(ard_initial_fitP,SpKf(j,:));   
    intafG=cumtrapz(ard_initial_fitP,SpGf(j,:));
    Ksatp(j,k)=Kh./(1+intafK(end));  
    Gsatp(j,k)=Gh./(1+intafG(end));  
    end   
%     plot(ard_initial_fitP,c_initialP);
%     figure
end

Kdrym=ones(nf,3)*Kh;
Gdrym=ones(nf,3)*Gh; 
Ksatp=[Ksatp Kdrym];  %%%K_mf  后面可以用Biot理论继续分析
Gsatp=[Gsatp Gdrym];

for i = 1:nP
Ksatp(125:end,i) = Ksatp(120,i);
Gsatp(125:end,i) = Gsatp(120,i);
end
save('Ksatp.mat','Ksatp');
save('Gsatp.mat','Gsatp');

