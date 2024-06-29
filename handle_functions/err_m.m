function err = err_m(v1,v2,varargin);
sv1 = size(v1); sv2 = size(v2);

if sv1 ~= sv2;
    error('The arrays have different dimensions');
end
Lv = length(sv1);

if nargin == 3
    if strcmp(varargin ,'tot')
        v1 = v1(:); v2 = v2(:); Lv = 1;
        err = max(2*abs(v1-v2)./(abs(v1)+abs(v2)),[],Lv);
    end
    if strcmp(varargin{1} ,'max')
        err = max(2*abs(v1-v2)./(abs(v1)+abs(v2)),[],Lv);
    end

elseif nargin == 4;
    if strcmp(varargin{1} ,'max')
        Lv = varargin{2};
        err = max(2*abs(v1-v2)./(abs(v1)+abs(v2)),[],Lv);
    end
else
    err = (2*abs(v1-v2)./(abs(v1)+abs(v2)));
end


end

