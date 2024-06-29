%local maxima of unidimensional real function .It can be used iteratively to find envelopes

function [val_m,varargout] = loc_max(x,varargin)

nx = length(x); dx  = diff(x); val_m =x;
if nargin > 1
    Nr = varargin{1};
else
    Nr = 1;
end
ind_tot = 1:nx;

for in = 1:Nr
    clear v;
    if in == 16;
        xxxx = 1;
    end
    nx = length(val_m);
    dx  = diff(val_m);
    x2 = val_m; k2 = 1;
    
    for k1 = 1: nx-2
        if(dx(k1)<0 && dx(k1+1)>0)
            v(k2) = k1+1;
            k2 = k2+1;
        end
    end
    
    if not(exist('v'));
        ix0 = find(dx == 0); Nx0 = length(ix0);
        indx = 1; v = zeros(Nx0,1); ivi = 1;
        while indx <= Nx0;
            ix = ix0(indx);
            if ix == 1
                v(ix) = 1; ivi = ivi+1;
            else
                if dx(ix-1) <0
                    ixf = ix;
                    while ismember(ixf+1,ix0); indx = indx+1; ixf = ixf+1;  end
                    
                    if in == 16 && indx == 123
                        xxx = 1;
                    end
                    
                    if (ixf == length(dx))
                        if(ixf > ix)
                            ivf = ivf - 1;
                        else
                            v(ivi) = ix;
                        end
                    else
                        
                        if (dx(ixf+1) > 0)
                            ivf = ivi + ixf-ix; v(ivi:ivf) = (ix:ixf); ivi = ivf+1;
                        end
                    end
                end
            end
            indx = indx+1;
        end
        v(ivi:Nx0) = [];
    end
    
    x2(v) = [];
    indv = 1:nx; indv(v) = []; val_m = x2;
    ind_tot = ind_tot(indv);
end

varargout{1} = ind_tot;

end

