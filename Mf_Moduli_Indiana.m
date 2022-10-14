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
por=0.139;   %
etaw=1e-3;
etagly=1e-0;
Kgly=4.36e9;
p=[3.6 5 10 20]'*1e6;
p_fit=[5 10 20]'*1e6;
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

Kdry=[7.647 14.7588 26.62 32.5876]'*1e9;
Kdry_fit=[14.7588 26.62 32.5876]'*1e9;
Kh=35.5e9;
Gh=26.6e9;
vh=(3*Kh-2*Gh)/(6*Kh+2*Gh);
Kex=Kh./Kdry_fit-1;
fex=fit(p_fit,Kex,'exp1');
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
Ksatp(130:end,i) = Ksatp(120,i);
Gsatp(130:end,i) = Gsatp(120,i);
end
save('Ksatp_india.mat','Ksatp');
save('Gsatp_india.mat','Gsatp');



