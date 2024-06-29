function Y = func_cell(X_cell,varargin); % It is linked to "precision_set2" function, and it calculates the double estimation from int
    X = X_cell{1}; sX = size(X); LX  = length(sX);
    DX = X_cell{3}; sDX = size(DX); LDX  = length(sDX);
    for ix = 1:LX; eval(['Ind_',num2str(ix),' = 1:sX(ix);']); end
    for ix = 1:LDX; eval(['Ind_D',num2str(ix),' = 1:sDX(ix);']); end
    num_arg2 = length(varargin);
    
    iarg = 1;
    while num_arg2 >= iarg;
        eval(['Ind_',num2str(iarg),' = varargin{',num2str(iarg),'};']); iarg = iarg+1;
    end
    
    iarg = 1;
    while num_arg2 >= iarg & iarg <= LDX;
        if not(eval(['isequal(Ind_D',num2str(iarg),',1)']));
        eval(['Ind_D',num2str(iarg),' = varargin{',num2str(iarg),'};']); end;  iarg = iarg+1;
    end
    
    str_indx = '(Ind_1'; for ix = 2:LX; str_indx = [str_indx,',Ind_',num2str(ix)]; end; str_indx = [str_indx, ')'];
    str_indD = '(Ind_D1'; for ix = 2:LDX; str_indD = [str_indD,',Ind_D',num2str(ix)]; end; str_indD = [str_indD, ')'];
    eval(['Y = (X_cell{2}',str_indD,' + double(X_cell{1}',str_indx,').*X_cell{3}',str_indD,'./X_cell{4});']);
end




