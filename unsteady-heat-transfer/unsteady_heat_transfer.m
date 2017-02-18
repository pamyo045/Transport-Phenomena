function y = unsteady_heat_transfer

PI    = pi();
% alfa  = 1;    %
% t     = 1;    % 
% b     = 1;    % m; centerSlab - surfSlab = halfWidthSlab
atb2  = 0.01;
% x     = linspace(0,b,10);    %
xb    = linspace(0,1,100);
% iTheta= zeros(1,length(xb));   %
nTheta= zeros(1,length(xb));    %
nT = zeros(1,length(xb));
% n     = 0;    %
eps   = 1e-20;

% atb2  = alfa*t/b^2
% xb    = x/b
% min(abs(curVal(abs(curVal)>0)));

% while min(abs(curVal(abs(curVal)>0))) > eps
%     curVal = (-1)^n/((n+0.5)*PI).*exp(-(n+0.5)^2*PI^2*alfa*t/b^2).*cos(n+0.5).*PI.*x./b;
%     theta  = theta + curVal;
%     n      = n + 1;
% end
% 
for n = 0:1e2
   oTheta    = nTheta;
   iTheta    = (-1)^n/((n+0.5)*PI).*exp(-(n+0.5)^2*PI^2*atb2).*cos(n+0.5).*PI.*xb;
   nTheta    = nTheta + iTheta;

%    if any(abs(oTheta - nTheta)) < eps
%        break;
%    end
end

nTheta  = 2*nTheta;
% figure; 
% hold on;
% for atb2 = [0.01 0.04 0.1 0.2 0.4 0.6 1.0];
%     for n = 1:1e4
%     %    oTheta    = nTheta;
%        Tleft  = 20;
%        Tright = Tleft;
%        Tini   = 100;
%        Cn     = 2./(n.*PI)*((-1).^n.*Tright - Tleft) + 2.*Tini.*(1 - cos(n.*PI))./(n.*PI);
%        iT     = Cn.*exp(-(n.*PI).^2.*atb2).*sin(n.*PI.*xb);
%        nT     = nT + iT;
% 
%     %    if any(abs(oTheta - nTheta)) < eps
%     %        break;
%     %    end
%     end
%     nT = Tleft + xb*(Tright - Tleft) + nT;
%     nT = (Tleft - nT)./(Tleft - Tini);
% %     nT = (nT - Tini)./(Tleft - Tini);
%     plot(2*xb,nT);
% end
% hold off;
% xlim([0 1]);
% n0      = 0:1e2-1;
% atb20   = 1;
% fh      = @(n)(-1).^n./((n+0.5).*PI).*exp(-(n+0.5).^2.*PI^2.*atb2).*cos(n+0.5).*PI.*atb20;
% figure;
% plot(n0,fh(n0));

% y.n     = n0;
% y.n     = y.n(1:10);
% y.theta = fh(n0);
% y.theta = y.theta(1:10); %(y.theta > 0.00001);
% y.one   = (-1).^(n0);
% y.one   = y.one(1:10); %(abs(y.one) > 0.00001);
% y.exp   = exp(-(n0+0.5).^2.*PI^2.*atb20);
% y.exp   = y.exp(1:10); %(abs(y.exp) > 0.00001);
% y.cos   = cos(n0+0.5);
% y.cos   = y.cos(1:10); %(abs(y.cos) > 0.00001);
% y.sum   = sum(y.theta);
% y.ntetha= nTheta(end);

% sum(fh(0:1e2))

figure;
plot(xb,nTheta);

% xlim([0 1]);
% ylim([0 1]);

end