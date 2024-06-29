function [x1 , no_itr, norm1,LAMBDA] = grad_desc(x0,F,DF,no_itr,error)
% nargin = no. of input arguments
if nargin <5 , no_itr = 20 ; end
if nargin <4 , error = 10^-5;no_itr = 20 ; end
if nargin <3 ,no_itr = 20;error = 10^-5; x0 = [1;1;1]; end

x1 = x0; loss_test = 0; Nk2 = length(x0); Nk = round(Nk2/2);
Fx1 = F(x1); Fold = Fx1; xold = x1; xold2 = x1; Fold2 = Fx1;
i = 0; lambda = 1e-3; LAMBDA(1) = lambda; change = 1; rate_impr = sqrt(7)/2; i2 = 1;
it_ret =1;
while i <= no_itr
    DFx1 = DF(x1);
    x1 = x1 - lambda.*DFx1; %x1 = x1.*(1 + 3e-2*(2*rand(Nk2,1)-1));
    %         r1 = x1(1:Nk);
    %         r1(r1 < 0) = 0; x1(1:Nk) = r1;
    Fx1 = F(x1);

 i == 500;
    if Fx1 < Fold
        lambda = lambda*rate_impr;
    else
        lambda = lambda./((rate_impr).^20)./1.1; x1 = xold; Fx1 = Fold;  xold = xold2; Fold = Fold2; 
    end


    i2 = i2+1; LAMBDA(i2) = lambda;
    xold2 = xold; Fold2 = Fold;
    xold = x1; Fold = Fx1;

    i = i + 1 ;
    norm1(i) = norm(Fx1);
    %     disp(['Iter ', num2str(i), ': loss = ', num2str(norm1(i))]);
    %         if norm1(i) < error , return , end
    if err_n(loss_test,Fx1) < 1e-5
        iterx = iterx+1; 
    else
        iterx = 0;
    end
    
    loss_test = Fx1;
    
    if i > no_itr || norm1(i) < error || iterx == 500 
        return ; 
    end

end
if 1  
xxx = 1;
end

end