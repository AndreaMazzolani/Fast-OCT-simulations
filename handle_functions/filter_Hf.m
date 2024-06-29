function Fx = filter_Hf(z,h)
% Needs to be cmopleted... fx and fxz needs to be the same

Nxh = length(h); Nxh2 = round(Nxh/2); tx = [1:Nxh] - Nxh2;  
Nyh = length(h); Nyh2 = round(Nyh/2); ty = [1:Nyh] - Nyh2;  
[Nxz,Nyz] = size(z); Nxz2 = round(Nxz/2);  Nyz2 = round(Nyz/2); 
txz = [1:Nxz] - Nxz2; tyz = [1:Nyz] - Nyz2;
[H_tmp,fx,fy] = gen_fft2D(+1,tx,ty,h); vecx = [fx(1), fx(end)]; vecy = [fy(1), fy(end)];
[H_tmp,fx,fy] = gen_fft2D(+1,tx,ty,h,vecx,vecy,Nxz,Nyz);
[Z_tmp,fxz,fyz] = gen_fft2D(+1,txz,tyz,z);

Fx = gen_fft2D(-1,fxz,fyz,H_tmp.*Z_tmp,txz,tyz);

end