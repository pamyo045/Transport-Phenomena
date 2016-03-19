function out = USHT2(tend,alfa,width)
% temporal domain
nt   = 2;
tau  = [0 tend*alfa/width^2];
dtau = abs(tau(2) - tau(1));

% spatial domain
nz   = 200;
eta  = linspace(-1,1,nz);
deta = abs(eta(2) - eta(1));

r    = dtau/deta^2;

% initial conditions
theta0 = 1; % @tau = 0

% boundary conditions
g1  = 0;    % @eta = -1;
g2  = 0;    % @eta = +1;

% initialize
theta  = zeros(nt,nz);
A = zeros(nz,1);
B = zeros(nz,1);
C = zeros(nz,1);
D = zeros(nz,1);

for j = 1:nt
    if(j == 1)
        theta(j,:)  = theta0;
        thetaOld    = theta0*ones(nz,1);
    else
        thetaOld   = theta(j-1,:)';
    end
    
    for i = 1:nz
        if(i == 1)      % lower boundary            
            A(i) =  0;
            B(i) =  1;
            C(i) =  0;
            D(i) =  g1;         
        elseif(i == nz) % upper boundary            
            A(i) =  0;
            B(i) =  1;
            C(i) =  0;
            D(i) =  g2;            
        else            % interior nodes  
            A(i) = -r;
            B(i) =  1 + 2*r;
            C(i) = -r;
            D(i) =  thetaOld(i);
        end  
    end 
        
    thetaNew   = TDMA(nz,A,B,C,D);
    theta(j,:) = (g1 - thetaNew)./(g1 - theta0);
end  
theta = 1 - theta(end,:);

out.Y = theta;
out.X = eta;
figure;
plot(eta, theta);
xlim([0 1]);
ylim([0 1]);
xlabel('\eta');
ylabel('\theta(t,z)');

% figure;
% x = eta.*width;
% % t = tau.*width^2./alfa;
% plot(x, theta);
% xlim([0 width]);
% % ylim([0 1]);
% xlabel('x (m)');
% ylabel('\theta(t,z)');
end

function E = TDMA(L,A,B,C,D)
    E     = zeros(L,1);
    BETA  = zeros(L,1);
    GAMMA = zeros(L,1);
    
    BETA(1)  = B(1);
    GAMMA(1) = D(1)/BETA(1);
    
 
   for I=2:L
        BETA(I)  = B(I) - A(I)*C(I - 1)/BETA(I - 1);
        GAMMA(I) = (D(I) - A(I)*GAMMA(I - 1))/BETA(I);
    end

    E(L) = GAMMA(L);
    for K = 1:L-1
        I    = L - K;
        E(I) = GAMMA(I) - C(I)*E(I + 1)/BETA(I);
    end
end