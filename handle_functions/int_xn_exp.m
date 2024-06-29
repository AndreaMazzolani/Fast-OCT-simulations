% definite integral of [  x^n exp(-(ax+b))^2   ]
% case in which the AN is underfilled.

function y = int_xn_exp(n,a,b);

if (a == 0)
    error('error, a = 0');
    else

    En = floor(n); En2 = floor(n/2);
    Eb1 = 0;

    for k = 0:En2
    Eb1 = Eb1 + factorial(n).*((-b).^(n-2.*k))./((4.^k).*factorial(k).*factorial(n-2*k));
    end

       y = sqrt(pi)./2./(a.^(n+1)).*Eb1;

    end
end
