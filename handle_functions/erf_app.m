function erf = erf_app(x)

if abs(x) < 1e-13
    erf = 0;
else

    tpi = 1/tan(pi/6); absx = abs(x); rex = real(x); re_imx = abs(real(x)./imag(x));
    condx1 = (absx <= 4); condx2 = (rex >= 0)&(absx>4)&(re_imx >= tpi); condx3 = (rex<0)&(absx>4)&(re_imx>=tpi); condx4 = (absx>=4)&(re_imx<tpi); anx1 = x(condx1);  anx4 = x(condx4);

    rx = real(anx1); ix = imag(anx1);
    ptmp = 0.3275911; ab1 = 0.254829592; ab2 = -0.284496736; ab3 = 1.421413741; ab4 = -1.453152027; ab5 = 1.061405429;
    coef = [ab5,ab4,ab3,ab2,ab1,0]; tx = 1./(1+ptmp.*abs(rx));
    p1 = 0; for kj = 0:5; p1 = p1+ coef(5-kj+1).*(tx).^kj; end

    ff1 =  ((1-p1.*exp(-rx.^2)).*sign(rx)) + exp(-rx.^2)./2./pi./rx.*(1-cos(2.*rx.*ix)+1i.*sin(2.*rx.*ix));
    for kk = 1:20
        ff1 = ff1 + 2/pi.*exp(-rx.^2).*exp(-kk^2/4)./(kk^2+4.*rx.^2).*(2.*rx.*(1-cos(2.*rx.*ix).*cosh(kk.*ix))+kk.*sin(2.*rx.*ix).*sinh(kk.*ix)+1i.*(2.*rx.*sin(2.*rx.*ix).*cosh(kk.*ix)+kk.*cos(2.*rx.*ix).*sinh(kk.*ix)));
    end
    erf = zeros(size(x));
    erf(condx1) = ff1;
    erf(condx2) = 1;
    erf(condx3) = (-1);
    erf(condx4)= ( sign(real(anx4))+ +exp(1).^((-1).*(1./anx4).^(-2)).*((-1).*pi.^(-1/2).*(1./anx4)+(1/2).*pi.^(-1/2).* ...
        (1./anx4).^3+(-3/4).*pi.^(-1/2).*(1./anx4).^5+(15/8).*pi.^(-1/2).*(1./anx4).^7+(-105/16).*pi.^( ...
        -1/2).*(1./anx4).^9+(945/32).*pi.^(-1/2).*(1./anx4).^11+(-10395/64).*pi.^(-1/2).*(1./anx4).^13));
end

end