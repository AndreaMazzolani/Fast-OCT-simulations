function err = err_n(v1,v2);
sv1 = size(v1); sv2 = size(v2);

Lv = length(sv1);
if (sv1 ~= sv2);
    error('The arrays have different dimensions');
end

v1 = v1(:); v2 = v2(:);
err = 2*norm(v1-v2)./(norm(v1)+norm(v2));
end