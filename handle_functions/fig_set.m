function [m_low,m_sup] = fig_set(B, varargin)

if nargin == 1 
    ns = 0.05; 
elseif nargin == 2
    ns = varargin{1};
else 
    error('Too many input arguments!');
end

[h,edges] = histcounts(B);
ch = cumsum(h);
Lh = length(h);

s=0; i=1;

i_low = max(find(ch <= ns*ch(Lh)));
i_sup = min(find(ch >= (1-ns)*ch(Lh)));

if isempty(i_low)
    i_low = 1;
end
if isempty(i_sup)
    i_sup = Lh;
end

m_low = edges(i_low);
m_sup = edges(i_sup+1);

end
