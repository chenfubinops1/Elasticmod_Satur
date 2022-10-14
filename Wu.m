function [p,q]=Wu(Ks,Gs,Kf,Gf,ar)    
            indx=find(ar==1); ar(indx)=0.99*ones(size(indx));

            obdx=find(ar<1);
            theta(obdx)=(ar(obdx)./((1-ar(obdx).^2).^(3/2))).*...
             (acos(ar(obdx)) -ar(obdx).*sqrt(1-ar(obdx).^2));
            fn(obdx)=(ar(obdx).^2./(1-ar(obdx).^2)).*(3.*theta(obdx) -2);

            ns=(3*Ks-2*Gs)./(2*(3*Ks+Gs));
            a=Gf./Gs -1; 
            b=(1/3)*(Kf./Ks -Gf./Gs); 
            r=(1-2*ns)./(2*(1-ns)); 
%             theta=(ar./((1-ar.^2).^(3/2))).*...
%              (acos(ar) -ar.*sqrt(1-ar.^2));
%             fn=(ar.^2./(1-ar.^2)).*(3*theta-2);          
            
            f1=1+a.*((3/2).*(fn+theta)-r.*((3/2).*fn+(5/2).*theta-(4/3)));
            f2=1+a.*(1+(3/2).*(fn+theta)-(r/2).*(3.*fn+5.*theta))+b.*(3-4*r);
            f2=f2+(a/2).*(a+3.*b).*(3-4.*r).*(fn+theta-r.*(fn-theta+2.*theta.^2));

            f3=1+a.*(1-(fn+(3/2).*theta)+r.*(fn+theta));

            f4=1+(a./4).*(fn+3.*theta-r.*(fn-theta));

            f5=a.*(-fn+r.*(fn+theta-(4/3))) + b.*theta.*(3-4*r);

            f6=1+a.*(1+fn-r.*(fn+theta))+b.*(1-theta).*(3-4.*r);

            f7=2+(a./4).*(3.*fn+9.*theta-r.*(3.*fn+5.*theta)) + b.*theta.*(3-4.*r);

            f8=a.*(1-2.*r+(fn./2).*(r-1)+(theta./2).*(5.*r-3))+b.*(1-theta).*(3-4.*r);

            f9=a.*((r-1).*fn-r.*theta) + b.*theta.*(3-4.*r);

            p=3*f1./f2; 
            q=(2./f3) + (1./f4) +((f4.*f5 + f6.*f7 - f8.*f9)./(f2.*f4));

            p=p./3; 
            q=q./5;  
end