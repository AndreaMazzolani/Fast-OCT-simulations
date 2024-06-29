function [x_prec,min_x,D_x,D_u] = precision_set(x,int_prec,varargin)

vec_prec = [8,16,32,64]; int_str = '';
for it = 1:3; if int_prec == vec_prec(it); int_str = num2str(int_prec); end; end

sx = size(x);
if isempty(int_str); error('The second input must be an integer among "8", "16", "32"'); end
str_int = ['uint',int_str]; sub_ind = 1:length(sx); num_var = length(varargin); 
if num_var == 1; sub_ind = varargin{1}; end


min_x = min(x,[],sub_ind);
max_x = max(x,[],sub_ind);
D_x = max_x-min_x;  

D_u = 2^int_prec - 1;

eval(['x_prec = ',str_int,'(D_u.*(x-min_x)./D_x);']); x_prec(D_x == 0) = 0; 


end


% TEST
% A= rand(1000,500,30); whos('A');

% [x_prec,x_min,D_x,D_u] = precision_set(A,8,[2,3]);


% A2 = x_min + double(x_prec).*D_x./D_u;
% err_n(A,A2)


