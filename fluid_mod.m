function Mif = pore_compress( BMf, SMf, ar)
%This function calculates effective compression modulus of the
% thin circular interlayer according to Tsai&Lee, 1998, IJSS.
% Basic assumptions:
% 1. Bondaries of the interlayer stays parallel during compression
% 2. Profile of the vertical lines becomes parabolic
% BMf, SMf = pore in-filler
% ar = aspect ratio
LMf = BMf-2/3*SMf;
alfaOB = sqrt(12*SMf./(LMf+2*SMf))./(2*ar);



Mif_temp=...
    LMf+2*SMf-(LMf.^2)./((LMf+2*SMf).*(alfaOB.*besseli(0,alfaOB))./...
    (2*besseli(1,alfaOB))-SMf);
% Because of Matlab limitations in handling high order values we have to
% rewrite the function and use a limiting value for the Bessel's functions 
% ratio and freqeuncy --> Inf
if ~isnan(Mif_temp)
    Mif = Mif_temp;
else
    Mif = LMf+2*SMf-(LMf.^2)./((LMf+2*SMf).*(1+1/8*alfaOB.^2)-SMf);
%     Mif = LMf+2*SMf-1/2*(LMf.^2)./((LMf+2*SMf).*alfaOB-SMf);
end