function A_noise = noise_add(A,varargin)
sA = size(A);

str_numb = 'R'; indmax = 1;
num_arg= length(varargin); str_noise = 'rand'; perc_noise = 1e-2;
if num_arg >= 1; perc_noise = varargin{1}; end
if num_arg >= 2; str_numb = varargin{2}; end
if num_arg >= 3; str_noise = varargin{3}; end
if num_arg >= 4; indmax = varargin{4}; end
if num_arg >= 5; error('Too many input data'); end

if strcmp(str_numb,'C'); randA = 2.*(rand(sA)+1i.*rand(sA))-1; 
elseif strcmp(str_numb,'R'); randA = 2.*(rand(sA))-1; 
    else; error('The second input argument must be "R" or "C" (real or complex)!');
end

if strcmp(str_noise,'rand');
    A_noise = A.*(1+perc_noise.*randA);
elseif strcmp(str_noise,'randmax');
    A_noise = A + perc_noise.*randA.*max(abs(A),[],indmax);
else
    error('Need to add other noises');
end

end
