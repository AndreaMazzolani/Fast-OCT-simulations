% This function calculates the displacement of the scatterers in z-direction after the mechanical loading, in presence of series of trains. 

%INPUT:
    % "alpha": Dimension (Nalpha,1). Array of strains. each strain defines a layer. "Nalpha" is hte number of strains.
    % "borders": Dimension (Nalpha +1,1). Array of all the borders of the layers of strains  alpha1 --> [b1,b2], alpha2 --> [b2,b3], alpha_j --> [b_j,b_(j+1)]
    % "z_U" Dimension (Nv,1). Array of the "unloaded" z-coordinates  
    % "varargin{1} = initial_disp". It is the initial bulk displacement of the scatterers. If it is missed then is set as "0", which means that the shallowest border is fixed. This displacement does not depend on the array of strains.

%OUTPUT:
    % "z_L" New positon of the scatterers z-coordinates after the meccanical loading
    % "varargout{1} = Disp". Dimension(Nalpha+1,1): Array of the borders after the mechanical loading. 


function [z_L, varargout] = load_scat_disp(alpha,borders,z_U,varargin)


LB = length(borders);
LS = length(alpha);

if not(issorted(borders))
    error(['The array of borders must be sorted!'])
end

if ((borders(1) > min(z_U)) | (borders(LB) < max(z_U)))
    error("The borders interval must contain the range of scatterers: bord(1) <= z_scat <= bord(end)");
end

Nz = length(z_U);
Disp = zeros(LB,1);
z_L = zeros(Nz,1);

z_bord = linspace(borders(1),borders(LB),LB*100);
alpha_bord = zeros(LB*100,1);

if nargin == 4
    initial_disp = varargin{1};
else
    initial_disp = 0;
end

Disp(1) = borders(1) + initial_disp;

if LB == LS+1
    for i_b = 2:(LB)
        izL = find((z_U < borders(i_b)) & (z_U >= borders(i_b-1))); % scatterers indices related to the i-th layer [borders(i_b-1),borders(i_b)]
        z_L(izL) = Disp(i_b-1) + (z_U(izL)-borders(i_b-1)).*(1-alpha(i_b-1)); % new scatterer positions of the i-th layer after the loading
        Disp(i_b) = Disp(i_b-1) + (borders(i_b)-borders(i_b-1)).*(1-alpha(i_b-1)); % Displacement of the borders after the loading
        iza = find((z_bord <= borders(i_b)) & (z_bord > borders(i_b-1))); % "z-bord" indices related to the i-th layer [borders(i_b-1),borders(i_b)]
        alpha_bord(iza) = alpha(i_b-1);
    end
else
    error(["The number of borders must equal the number of the strains + 1  -->  (|Borders| = |Strains| +1)"]);
end

varargout{1} = Disp; 
varargout{2} = z_bord;
varargout{3} = alpha_bord; 

end
