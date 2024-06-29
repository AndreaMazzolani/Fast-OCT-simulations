function output = multi_indices_composition(X,Ij,j)

sx = size(X); Nx = length(sx);
sj = size(Ij);

if j == 1; sjx = sx(2:Nx);
elseif j == Nx; sjx = sx(1:(Nx-1));
else; sjx = [sx(1:(j-1)) ,sx((j+1):Nx)];
end
if not(isequal(sj,sjx)); error('The dimensions do not match!'); end;

str_output = []; str_input = [];

for ix = 1:Nx
    if not (ix == j)
        nix = num2str(ix);
        str_output = [str_output,'R',nix,','];
        str_input = [str_input,'1:sx(',nix,'),'];
    end
    if ix < j;  str_out_pre = str_output; end
end

Nx_out_pre = length(str_out_pre); Nx_out = length(str_output); Nx_in = length(str_input);
str_out_post = str_output(Nx_out_pre+1:Nx_out-1);

Nx_in = length(str_input);
str_output = str_output(1:Nx_out-1);
str_input = str_input(1:Nx_in-1);


str_out_pre = str_out_pre(1:Nx_out_pre-1);
% tic;
eval(['[',str_output,'] = ndgrid(',str_input,');']);

% toc;

if j == 1; eval(['output = X(sub2ind(sx,Ij,',str_output,'));']);
elseif j == Nx; eval(['output = X(sub2ind(sx,',str_output,',Ij));']);
else eval(['output = X(sub2ind(sx,',str_out_pre,',Ij,',str_out_post,'));']); end

% [rows, columns] = ndgrid(1:size(Y, 1), 1:size(Y, 2));
% output = X(sub2ind(size(X), Y, rows, columns))


Test = 0 ; % Test to run in another script to test the function;
if Test


    N1 = 5; N2 = 3; N3 = 7; N4 = 5; output1 = zeros(N1,N2,N4);
    X = rand(N1,N2,N3,N4); j = 3; Ij = ceil(rand(N1,N2,N4).*N3); % Ij is an array of indices for the j-th dimension, for each other dimension!

    for m1 = 1:N1
        for m2 = 1:N2
            for m4 = 1:N4
                output1(m1,m2,m4) = X(m1,m2,Ij(m1,m2,m4),m4);
            end
        end
    end

    output2 = multi_indices_composition(X,Ij,j);
    max(abs(output1(:) - output2(:)))
end


end



