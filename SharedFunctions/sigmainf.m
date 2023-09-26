function [SUM]=sigmainf(arf,k,rho,phi,theta,mode,limit)
%SIGMAINF   用于求解正负无穷求和的收敛值
%   k 波矢 value
%   rho & phi 柱坐标 mat
%   theta 微元 value
%   mode 选择求和公式的模型 可选"1，2，3，4，5"
%   limit 选择的上下限 任意整数 如"100, 200, 1000"
%

%   Copyright 2023/09/15 Lusen.

if mode == "1"
    Sum = 0;
    for n = -limit:limit
        temp=(arf*cos(arf*pi)+1i*n*sin(arf*pi))/(arf^2-n^2).*(1i)^(n+1).*exp(1i.*n.*phi).*sin(arf.*pi).*besselj(n,k.*rho.*sin(theta)); % A1
        Sum = Sum + temp;
    end
end

if mode == "2"
    Sum = 0;
    for n = -limit:limit
        temp=((arf-2)*cos(arf*pi)+1i*n*sin(arf*pi))/((arf-2)^2-n^2).*(1i)^(n+1).*exp(1i.*n.*phi).*sin(arf.*pi).*besselj(n,k.*rho.*sin(theta)); %A2
        Sum = Sum + temp;
    end
end

if mode == "3"
    Sum = 0;
    for n = -limit:limit
        temp=(n*cos(arf*pi)+1i*arf*sin(arf*pi))/(arf^2-n^2).*(1i)^(n).*exp(1i.*n.*phi).*sin(arf.*pi).*besselj(n,k.*rho.*sin(theta)); % B1
        Sum = Sum + temp;
    end
end

if mode == "4"
    Sum = 0;
    for n = -limit:limit
        temp=(n*cos(arf*pi)+1i*(arf-2)*sin(arf*pi))/(n^2-(arf-2)^2).*(1i)^(n).*exp(1i.*n.*phi).*sin(arf.*pi).*besselj(n,k.*rho.*sin(theta)); %B2
        Sum = Sum + temp;
    end
end

if mode == "5"
    Sum = 0;
    for n = -limit:limit
        temp=((arf-1)*cos(arf*pi)+1i*n*sin(arf*pi))/((arf-1)^2-n^2).*(1i)^(n+1).*exp(1i.*n.*phi).*sin(arf.*pi).*besselj(n,k.*rho.*sin(theta)); %C1
        Sum = Sum + temp;
    end
end

SUM = Sum;


end

