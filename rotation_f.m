function [theta,x1,y1,x2,y2]=rotation_f(A,k_noise) % Input cropped matrix image A, output rotation angle, x-direction pixel difference, y-direction pixel difference
d=2;                  % The radius of a single light spot
k=0.0;               

f1=max(max(A));        
A=A/f1;             

Ad=(A>k);           
A=A.*Ad;
Ap=A';     
[f1b,f1a]=find(Ap==1,1);
B=A; 

  for i=-d:d
    for j=-d:d
      A(f1a+i,f1b+j)=0;
    end
  end                   % Spot area clear to 0

f2=max(max(A));
Ap=A';                 
[f2b,f2a]=find(Ap==f2,1);

  for i=-d:d
    for j=-d:d
      A(f2a+i,f2b+j)=0;
    end
  end                   
A2=A>0;
f3=max(max(A));
noise=sum(A(:))/sum(A2(:))*k_noise;
B=B-noise;
B(B<0)=0;

B1=B(f1a-d:f1a+d,f1b-d:f1b+d); 

B2=B(f2a-d:f2a+d,f2b-d:f2b+d);  

x1 = sum(sum(B1 .* [1 : size(B1, 2)])) /(sum(B1(:))+eps)+f1b-d;  
y1 = sum(sum(B1' .* [1 : size(B1, 1)])) /(sum(B1(:))+eps)+f1a-d; 
x2 = sum(sum(B2 .* [1 : size(B2, 2)])) /(sum(B2(:))+eps)+f2b-d;  
y2 = sum(sum(B2' .* [1 : size(B2, 1)])) /(sum(B2(:))+eps)+f2a-d;
t_theta= (x2-x1)/(y2-y1+eps); 
theta=atan(t_theta)*180/pi; % rotation angle
 if f3>(1+f2)*0.44
     theta=NaN;
 end                         %confidence level
end