function [x_prec_cell,varargout] = precision_set2(x,int_prec,varargin)

vec_prec = [8,16,32,64]; int_str = '';
for it = 1:3; if int_prec == vec_prec(it); int_str = num2str(int_prec); end; end

sx = size(x); Lx  = length(sx);
if isempty(int_str); error('The second input must be an integer among "8", "16", "32"'); end
str_int = ['uint',int_str]; sub_ind = 1:length(sx); num_var = length(varargin);
if num_var == 1; sub_ind = varargin{1}; end
if num_var == 2; error('Too many input'); end

min_x = min(x,[],sub_ind);
max_x = max(x,[],sub_ind);
D_x = max_x-min_x; DX = D_x.*ones(sx);

D_u = 2^int_prec - 1;

eval(['x_prec = ',str_int,'(D_u.*(x-min_x)./D_x);']); x_prec(DX == 0) = 0;

x_prec_cell = cell(4,1);
x_prec_cell{1} = x_prec;
x_prec_cell{2} = min_x;
x_prec_cell{3} = D_x;
x_prec_cell{4} = D_u;


f = @func_cell;


if nargout == 2; varargout{1} = f; end
end




% TEST:

% A = rand(30,151,117,48);  
 
%%%% To restore the array "A", call the function "func_cell"

% A_prec_cell = precision_set2(A,8,[2,3]);

% A_est = func_cell(A_prec_cell,3:7,20:40,1:5,13:17);
% err_n(A(3:7,20:40,1:5,13:17),A_est)

