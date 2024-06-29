
function [coef] = poly_lagrange(xnodi,ynodi)
xnodi = xnodi(:); ynodi = ynodi(:);
n = length(xnodi);
A = zeros(n);  znodi = zeros(n,1);

for i = 1:n
    A(:,i) = xnodi.^(i-1);
    znodi(i,:,:,:) = ynodi(i,:,:,:);
end

coef  = A\znodi;

end