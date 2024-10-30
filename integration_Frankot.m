function z = integration_Frankot(Ip,Iq,dx,dy)

z = integration_Frankot_mathias(Ip,Iq,dx,dy);

function Iz = integration_Frankot_mathias(Ip,Iq,dx,dy)

[lig,col] = size(Ip);


for i1=1:lig
    for j1=1:col
        u(i1,j1) = (j1-1)-col./2;
        v(i1,j1) = (i1-1)-lig./2;
    end
end
u = u./col.*(2*pi)./dx;
v = v./lig.*(2*pi)./dy;
u = fftshift(u);
v = fftshift(v);

tfp = fft2(Ip);
tfq = fft2(Iq);


%% integration
A  = -i.*( tfp.*u + tfq.*v );
B  = (u.^2 + v.^2);

B(find(B==0)) = eps;%1;%0.0000001; % ne pas diviser par 0
tmp = A ./ B;

tmp2 = real(ifft2(tmp));

[lig,col] = size(tfq);
tmpp = dx*real(tfp(1,1))/(lig*col);
tmpq = dy*real(tfq(1,1))/(lig*col);
for i1=1:lig
    for j1=1:col
        pentes(i1,j1) = (j1-1)*tmpp + (i1-1)*tmpq;
    end
end

Iz = pentes + tmp2;

function z=morel(p,q,varargin)       

[nbl,nbc]=size(p);
dx=1;
dy=1;

if length(varargin)>= 1
   dx=varargin{1};
   dy=dx;
end
if length (varargin) == 2
   dy=varargin{2};
end

tf_p=fft2(p);
tf_q=fft2(q);

u=fftshift(((0:nbc-1)-nbc/2)/nbc)/dx*2*pi;
v=fftshift(((0:nbl-1)-nbl/2)/nbl)/dy*2*pi;
for k=1:nbl
   U(k,:)=u;
end
for k=1:nbc
   V(:,k)=v';
end


tfZ=-i*(U.*tf_p+V.*tf_q)./(U.^2+V.^2+eps);

pente_x=real(tf_p(1,1))/(nbc*nbl);
pente_y=real(tf_q(1,1))/(nbc*nbl);

x=0:dx:dx*(nbc-1);
y=0:dy:dy*(nbl-1);

for l=1:nbl
   for c=1:nbc
      const(l,c)=pente_x*x(c)+pente_y*y(l);
   end
end
z=real(ifft2(tfZ))+const;
