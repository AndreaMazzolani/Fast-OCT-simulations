%fourier_ transform with Nt and Nk chosen.

function [g] = ff_tr(F,It,freq,refind);

Nt = length(It);
DT = zeros(Nt,1);
DT(1) = It(2)-It(1);
for it = 2:(Nt-1);
    DT(it) = (It(it+1)-It(it-1))/2;%DT(it-1)/2 + (It(it+1)-It(it))/2; ? I do not remember why I put that before
end
DT(Nt) = It(Nt)-It(Nt-1);

sx = size(F);
Lx = length(squeeze(sx));

if (Lx == 1 || sx(2) == 1)
    g = zeros(length(freq),1);
    for k = 1:length(freq);
        if length(It)< 3
            for j = 1: length(It);
                dt = DT(j);
                g(k) = g(k) + F(j)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;

            end
        else
            for j = 2: length(It)-1;
                dt = DT(j);
                g(k) = g(k) + F(j)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;
            end
            g(k) = g(k) + F(1)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(1)),2*pi))*DT(1)/2 + F(Nt)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(Nt)),2*pi))*DT(Nt)/2;
        end
    end
    g = g(:);
elseif (Lx == 2)
    g = zeros(length(freq),sx(2));
    for k = 1:length(freq);
        if length(It)< 3
            for j = 1: length(It);
                dt = DT(j);
                g(k,:) = g(k,:) + F(j,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;
            end
        else
            for j = 2: length(It)-1;
                dt = DT(j);
                g(k,:) = g(k,:) + F(j,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;
            end
            g(k,:) = g(k,:) + F(1,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(1)),2*pi))*DT(1)/2 + F(Nt,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(Nt)),2*pi))*DT(Nt)/2;
        end
    end
elseif (Lx ==3)
    g = zeros(length(freq),sx(2),sx(3));
    for k = 1:length(freq);
        if length(It)< 3
            for j = 1: length(It);
                dt = DT(j);
                g(k,:,:) = g(k,:,:) + F(j,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;
            end
        else
            for j = 2: length(It)-1;
                dt = DT(j);
                g(k,:,:) = g(k,:,:) + F(j,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;
            end
            g(k,:,:) = g(k,:,:) + F(1,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(1)),2*pi))*DT(1)/2 + F(Nt,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(Nt)),2*pi))*DT(Nt)/2;
        end
    end
elseif (Lx ==4)
    g = zeros(length(freq),sx(2),sx(3),sx(4));
    for k = 1:length(freq);
        if length(It)< 3
            for j = 1: length(It);
                dt = DT(j);
                g(k,:,:,:) = g(k,:,:,:) + F(j,:,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;
            end
        else
            for j = 2: length(It)-1;
                dt = DT(j);
                g(k,:,:,:) = g(k,:,:,:) + F(j,:,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(j)),2*pi))*dt;
            end
            g(k,:,:,:) = g(k,:,:,:) + F(1,:,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(1)),2*pi))*DT(1)/2 + F(Nt,:,:,:)*exp(-1i.*mod(2*pi*freq(k)*refind*(It(Nt)),2*pi))*DT(Nt)/2;
        end
    end
else
    error('Dimensioon too high')
end

xxx = 1;
end

