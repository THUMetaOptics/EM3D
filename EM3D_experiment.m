%% Image read in
close all; clear;
orial_pic = imread("raw_image.bmp");
pic1 = orial_pic(525:1824,981:2020);
%% Polarization calculation
[T,PHI,DoP] = T_PHI_DoP(pic1);
DoP = medfilt2(DoP,[5,5]);
PHI = medfilt2(PHI,[5,5]);                                       %denoising
PHI(PHI<0) = PHI(PHI<0)+2*pi;
figure(1);imagesc(pic1);title('Raw image'); colormap("gray")
figure(2);imagesc(DoP);title('DOLP');colorbar;caxis([0,0.2]);colormap(brewermap(168,'rdpu'));
PHI2 = PHI-pi/2;
PHI2(PHI2<0) = PHI2(PHI2<0)+2*pi;
PHI2(PHI2>pi) = PHI2(PHI2>pi)-pi;
figure(3);imagesc(PHI2);title('AOP');colorbar;caxis([0,pi]);colormap(brewermap(168,'ylorrd'));

%% Surface normal before correction
n=1.5;
elevation2 = abs(theta_cal_fun(DoP,n));
elevation2 = medfilt2(elevation2);
figure(4);imagesc(elevation2); stitle=sprintf('Raw zenith angle'); title(stitle);colorbar;colormap(brewermap(168,'Reds'));caxis([0,0.5*pi]);
PHI3 = PHI2+pi/2;
figure(5);imagesc(PHI3); stitle=sprintf('Raw azimuth angle'); title(stitle);colorbar;colormap(brewermap(168,'orrd'));caxis([0,2*pi]);
Gx = tan(elevation2).*cos(-PHI3);                             % x component for normal vector
Gy = tan(elevation2).*sin(-PHI3);                             % y component for normal vector
dx =0.025;
dy = 0.025;
z = integration_Frankot(Gx,Gy,dx,dy);

figure(6);
surfl(imrotate(z,90));
shading interp;
colormap gray(256);
colormap(brewermap(168,'*Reds'));
lighting phong;
grid off;
axis off;

%% Raw absolute depth from PSFs
[h,l] = size(orial_pic);
LA = 2000;
pa = 00;
LB = 1600;
pb = 290;
L = 400;
or = 4;
I2 = double(orial_pic(h/2-LA/2+pa+1:h/2+LA/2+pa,l/2-LB/2+pb+1:l/2+LB/2+pb));

I2 = medfilt2(I2,[5,5]);
for i = 1:L
    for j = 1:L
        
        H(i,j) = 1/4.*(1-cos(2.*pi.*i./(L-1))).*(1-cos(2.*pi.*j./(L-1)));
    end
end
H = imresize(H,[L,L]);                                               %Hanning window

rx = LA/L*or-or+1;
ry = LB/L*or-or+1;
angle = zeros(rx,ry);
anglec = zeros(rx,ry);
era = cell(rx,ry);
ima = cell(rx,ry);
for i = 1:(rx)
    for j = 1:(ry)                                       
        temp1 = I2((i-1)*L/or+1:(i+or-1)*L/or,(j-1)*L/or+1:(j+or-1)*L/or);    
        temp1 = imadjust(temp1/256);        
        mr = 1;
        Wh = temp1.*H;
        Wh = imresize(Wh,[mr*L,mr*L],'bilinear');
        ima{i,j} = temp1;
        Cep = fftshift(ifft2(log((abs(fftshift(fft2(Wh)))).^2)));      %cepstrum
        dp = 32*mr;
        pmin = 0.7.*dp/2;
        pmax = 1.3.*dp/2;
        w = 1;
        [m1,n1] = size(Cep);
        for m = 1:m1
            for n = 1:n1
                if sqrt((m-m1/2-1).^2+(n-n1/2-1).^2)/2 >= pmin&&sqrt((m-m1/2-1).^2+(n-n1/2-1).^2)/2 <= pmax
                    Cep(m,n) = Cep(m,n);
                else Cep(m,n) = 0;
                end
                if m > (m1/2-w+1) && m < (m1/2+w+1)
                    Cep(m,n) = 0;
                end
                if n > (n1/2-w+1)&&n < (n1/2+w+1)
                    Cep(m,n) = 0;
                end         
            end
        end                                                                 %cepstrum area of interest

        A = abs(Cep);
        era{i,j} = A;
        A = imgaussfilt(A,[2,2]);
         
        [t0,x1,y1,x2,y2] = rotation_f(A,1.5);                    %rotation angle                 
        if t0<0
            t0 = t0+180;
        end
        oar = abs(i-rx/1.6)/sqrt((rx/2).^2+(ry/2).^2)*3-abs(j-ry/1.5)/sqrt((rx/2).^2+(ry/2).^2)*7;         %Correction of optical imaging distoration, parameters obtained by camera calibration
        angle(i,j) = t0;
        anglec(i,j) = t0-oar;
    end
end
cc = nanmean(nanmean(anglec));
anglec(anglec>cc+12) = nan;
anglec(anglec<cc-12) = nan;
cc = nanmean(nanmean(anglec));
 anglec(anglec>cc+8) = nan;
 anglec(anglec<cc-8) = nan;                                    %eliminate outliers
DEDH = 9600./(-anglec+220);
figure(7),imagesc(DEDH);title('Raw PSF depth');set(imagesc(DEDH),'alphadata',~isnan(DEDH)*0.9+0.1);colormap(brewermap(168,'*GnBu'));caxis([57.5,62]); %Raw absolute depth

%% Completing and filtering absolute depth (PSF disparity)
 [fx,fy] = size(anglec);
anglecf = fillmissing(anglec,'movmean',45);
anglecf = fillmissing(anglecf.','movmean',6);
anglecf = anglecf.';
anglecf = imgaussfilt(anglecf,[4.7,1.3]);

%% Azimuth angle of completed absolute depth
depthgy = diff(anglecf);
depthgy = imresize(depthgy,[fx,fy]);
depthgx = diff(anglecf.').';
depthgx = imresize(depthgx,[fx,fy]);
for i = 1:fx
    for j = 1:fy 
        x = depthgx(i,j);
        y = depthgy(i,j);
        if x == 0
            if y > 0
                PHIDH(i,j) = 90;
                
            else PHIDH(i,j) = 270;               
            end        
        else
            if y >= 0
                if x > 0
                    PHIDH(i,j) = atand(y/x);
                else
                    PHIDH(i,j) = 180-atand(-y/x);
                end               
            else if x > 0
                    PHIDH(i,j) = 360-atand(-y/x);
                else PHIDH(i,j) = 180+atand(y/x);
                end
            end
        end
    end
end    
PHIDH = PHIDH*pi/180;
PHIDH(PHIDH < 0.6) = 2*pi;                                     %reduce 2*pi gap
figure(8),imagesc(PHIDH);caxis([0,2*pi]);title('PSF azimuth angle');colormap(brewermap(168,'orrd'));       %azimuth angle of completed absolute depth

%% Correction of polarization azimuth angle
RR = 1;
PHIDH = imresize(PHIDH,[LA/RR/2-350,LB/RR/2-280],'bilinear');
PHIC = PHI;
for i = 1:LA/RR/2-350
    for j = 1:LB/RR/2-280
        PHID(i,j) = abs(PHI(i,j)-PHIDH(i,j));
        if PHID(i,j) > pi/2&&PHID(i,j) < pi*3/2
         PHIC((i-1)*RR+1:i*RR,(j-1)*RR+1:j*RR) = PHI((i-1)*RR+1:i*RR,(j-1)*RR+1:j*RR)+pi;
        end
    end
end
PHIC(PHIC > 2*pi)=PHIC(PHIC > 2*pi)-2*pi;
Gxc = tan(elevation2).*cos(-PHIC);  
Gyc = tan(elevation2).*sin(-PHIC);
figure(9),imagesc(PHIC),caxis([0,2*pi]);title('Corrected polarization azimuth angle');colormap(brewermap(168,'orrd'));

%% Completed absolute depth and its zenith angle
DEDH = 9600./(-anglecf+220);                                %convert angles to depth, parameters obtained by camera calibration
AVD = mean(mean(DEDH));
PFOV = 3.520;                                                         %calibrated camera parameter
pw = AVD/PFOV/fx;
abdgy = diff(DEDH);
abdgy = imresize(abdgy,[fx,fy]);
abdgx = diff(DEDH.');
abdgx = imresize(abdgx.',[fx,fy]);
GABD = sqrt(abdgy.^2+abdgx.^2).*pw;
CM = zeros(fx,fy);
CM(anglec > 0) = 1;
GWC = GABD.*CM;
GWC(GWC == 0) = NaN;
GWC(GWC > 1.35) = 1.35;
figure(10),imagesc(DEDH);caxis([58,62]);title('Completed absolute depth');colormap(brewermap(168,'*GnBu'));         %Completed absolute depth
figure(11),imagesc(GWC);set(imagesc(GWC),'alphadata',~isnan(GWC)*0.92+0.08);title('PSF zenith angle');colormap(brewermap(168,'Reds'));caxis([0,1.57]);    %zenith angle of Completed absolute depth
CM2 = CM;
CM2(CM2 == 0) = nan;
AVD2 = nanmean(nanmean(DEDH.*CM2));

%% Correction of surface index
n = 1;
doplr = imresize(DoP,[fx,fy]);
for i = 1:10
    GDP = abs(theta_cal_fun(doplr,n));
    dif1(i) = abs(sum(sum((GDP-GABD).*CM)));
    if i>1
        if dif1(i) > dif1(i-1)
            n = n-0.1;
            break
        end
    end
    n = n+0.1;
end
n1 = n-0.1;
for i = 1:20
    GDP = abs(theta_cal_fun(doplr,n1));
    dif2(i) = abs(sum(sum((GDP-GABD).*CM)));
    if i>1
        if dif2(i) > dif2(i-1)
            n1 = n1-0.01;
            break
        end
    end
    n1 = n1+0.01;
end
n2 = n1-0.01;
for i = 1:20
    GDP = abs(theta_cal_fun(doplr,n2));
    dif2(i) = abs(sum(sum((GDP-GABD).*CM)));
    if i > 1
        if dif2(i) > dif2(i-1)
            n2 = n2-0.001;
            break
        end
    end
    n2 = n2+0.001;
end
%% Correction of polarization zenith angle
elevationc = abs(theta_cal_fun(DoP,n2));
Gxc = tan(elevationc).*cos(-PHIC);                           
Gyc = tan(elevationc).*sin(-PHIC); 
figure(12);imagesc(elevationc); colorbar;colormap(brewermap(168,'Reds'));caxis([0,1.57]);

 %% Depth after correction
dx = 0.025;
dy = 0.025;
z = integration_Frankot(Gxc,Gyc,dx,dy);                    %corrected polarization 3D

ABDM = (z-mean(mean(z)))/dx/PFOV/fx+AVD2;       %Refined 3D surface with absolute depth
figure(13);
surf(imrotate(ABDM,90));
shading interp;
colormap gray(256);
colormap(brewermap(168,'*GnBu'));
lighting phong;
grid off;
axis off;
stitle=sprintf('Surface recovery Frankot'); title(stitle);caxis([58,62]);

%% Depth accuarcy
zc = mean(ABDM(LB/4-50:LB/4+50,:));                     %section

px = pw;
LA = 1040;
y = -px*2.75+px*5.5/LA:px*5.5/LA*2:px*2.75;
pc = 263;
pr = 265;
for i = 1:520
    zg(i) = -sqrt(1-((i-pc)/pr).^2)*px*pr/90+62.16;      %ground truth
end   

meanerror = mean(abs(zg-zc))/mean(zg);                 %meanerror

figure(14),p1 = plot(y,zc,'--','color',[178/255,24/255,43/255],'LineWidth',1.3) ;hold on;
p2 = plot(y,zg,'LineWidth',1.3,'color',[33/255,102/255,172/255]);
xlabel('x space (cm)','FontSize',12,'Fontname','Arial'),ylabel('z space (cm)','FontSize',12,'Fontname','Arial'),xlim([-3.5,3.5]);ylim([58.4,63.2]);
legend([p1 p2], {'Calculated','Ground Truth'},'FontSize',20,'Fontname','Arial','Location','southeast')

 %%
 function [T,PHI,DoP] = T_PHI_DoP(orial_pic)
[m,n,k] = size(orial_pic);
odd_n = 1:2:n;
even_n = 2:2:n;
odd_m = 1:2:m;
even_m = 2:2:m;
pol_0 = orial_pic(odd_m,odd_n,1);
pol_45 = orial_pic(odd_m,even_n,1);
pol_135 = orial_pic(even_m,odd_n,1);
pol_90 = orial_pic(even_m,even_n,1);
I90 = im2double (pol_90);
I45 = im2double (pol_45);
I135 = im2double (pol_135);
I0 = im2double (pol_0);

[m,n,~] = size(I0);
T = (I0 + I90);
%DOLP
cc = (I90-I0).^2+(I45-I135).^2;
DoP = sqrt(cc)./(I0+I90);
DoP(T ==0) = 0;

%angle of polarization
PHI = (1/2)*atan((I45 - I135)./(I0-I90+0.00000001));


for i = 1:1:m
    for j = 1:1:n
        ind1 =  I90(i,j) < I0(i,j) ;
        ind2 = I45(i,j) < I0(i,j) ;

        if ind1
            PHI(i,j)=PHI(i,j)+pi/2;
        end
        if ind2
            PHI(i,j)=PHI(i,j)+pi;
        else
            PHI(i,j) = PHI(i,j);
        end
    end
end

[m1,n1,k1] = size(I0);
DoP=DoP-0.01;                                                     %correction of dark noise
[p_depth1,q_depth1] = gradient(I0);
p_afthercor1 = zeros(m1,n1);
q_afthercor1 = zeros(m1,n1);
end