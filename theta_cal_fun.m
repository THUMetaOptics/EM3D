function [theta1]=theta_cal_fun(P,n)

fenzi = 2*P +2*n^2.*P - 2*n^2+ n^4 +P.^2 + 4*n^2*P.^2-n^4*P.^2-4*n^3.*P.*sqrt(-(P-1).*(P+1)) +1;
fenmu = n^4*P.^2 + 2*n^4*P + n^4 +6*n^2*P.^2 +4*n^2*P - 2*n^2 + P.^2 + 2*P +1;
cos_theta1 = sqrt(fenzi./fenmu);
theta1 = acos(cos_theta1);

end
