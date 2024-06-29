function err = err_mn(v1,v2,varargin);
sv1 = size(v1); sv2 = size(v2);

if sv1 ~= sv2;
    error('The arrays have different dimensions');
end
%         err = 2*abs(v1-v2)./mean(abs(v1)+abs(v2),'all');

M1 = max(abs(v1),[],'all');
M2 = max(abs(v2),[],'all');
err = 2*abs(v1-v2)./min(M1,M2);

% Lv = length(sv1);
% if nargin == 3
%     if strcmp(varargin ,'tot')
%         v1 = v1(:); v2 = v2(:); Lv = 1;
%         err = 2*abs(v1-v2)./mean(abs(v1)+abs(v2),[],'all');
%     end
%     if strcmp(varargin{1} ,'max')
%         err = 2*abs(v1-v2)./mean(abs(v1)+abs(v2),[],Lv);
%     end
% 
% elseif nargin == 4;
%     if strcmp(varargin{1} ,'max')
%         Lv = varargin{2};
%         err = 2*abs(v1-v2)./mean(abs(v1)+abs(v2),[],Lv);
%     end
% else
%     err = (2*abs(v1-v2)./(abs(v1)+abs(v2)));
% end


end

