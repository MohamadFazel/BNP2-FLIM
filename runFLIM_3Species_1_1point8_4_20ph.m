
%Example script for simulated FLIM data similar to those in Fig. 2 in the 
%main text with lifetimes of 1ns, 1.8ns and 4.5ns
addpath('Functions')

%% Populating EmpParams structure
%EmpParam contain parameters of the experiment described in the following:
%These parameters are used in data generation and some of them are used in
%the algorithm.
%
% Xi_X:   Confocal X center (pixel)
% Yi_X:   Confocal Y center (pixel)
% Zi_X:   Confocal Z center (pixel)
% Mu:     Array of excitation rate of the species (1/ns)
% Tau:    Array of lifetimes of the species (ns)
% DeltaP: Pulse duration (ns)
% Sig_IRF:Sigma of IRF (ns)
% T_IRF:  Mean of IRF (ns)
% T:      Time interval between two consecurive excitation pulses (ns)
% OmegaX: HFWM of confocal PSF along X and Y axis (pixel)
% OmegaZ: HFWM of confocal PSF along Z axis (pixel)   
% PSF:    The combination of detection and illumination PSFs (inline function)   
%

M=6;
EmpParam.Xi_X = [];
XPixNum = 20;
YPixNum = 5;
EmpParam.PixelSize = 0.3922; %mu
for pp = 1:XPixNum
    EmpParam.Xi_X = cat(1,EmpParam.Xi_X,(pp-0.5)*ones([YPixNum,1]));
end
EmpParam.Xi_Y = repmat((0.5:YPixNum-0.5)',[XPixNum,1]);
EmpParam.Xi_Z = zeros(XPixNum*YPixNum,1);
EmpParam.Xi_X = EmpParam.PixelSize*EmpParam.Xi_X;
EmpParam.Xi_Y = EmpParam.PixelSize*EmpParam.Xi_Y;
EmpParam.Mu = 0.001*ones(3,1);
EmpParam.Tau = [1 1.8 4.5];
EmpParam.Dp = 0.1;
EmpParam.Sig_IRF = 0.8;
EmpParam.T_IRF = 12.21;
freq=79997407;
EmpParam.T = 10^9/freq;
EmpParam.OmegaX = 1.4*EmpParam.PixelSize; %mu
EmpParam.OmegaZ = 4*EmpParam.PixelSize; %mu
EmpParam.PSF = @(X,Y,Xp,Yp,OmegaX,OmegaZ) OmegaZ*exp(-2*((X-Xp).^2/OmegaX^2 + ...
    (Y-Yp).^2/OmegaX^2));

%%

RhoIn(1).Rho = @(X,Y,Z) 70*erfc((X-XPixNum*EmpParam.PixelSize/4));
RhoIn(2).Rho = @(X,Y,Z) 70*erfc((XPixNum*EmpParam.PixelSize-6.5*X/5));
RhoIn(3).Rho = @(X,Y,Z) 150*erfc(2*abs(X-XPixNum*EmpParam.PixelSize/2)/3);

[Xg,Yg] = meshgrid(EmpParam.PixelSize*(0.5:19.5),EmpParam.PixelSize*(0.5:5.5));
figure;surf(Xg,Yg,0.001*RhoIn(1).Rho(Xg,Yg,0))
hold;surf(Xg,Yg,0.001*RhoIn(2).Rho(Xg,Yg,0))
surf(Xg,Yg,0.001*RhoIn(3).Rho(Xg,Yg,0))

D = EmpParam.PixelSize; %Pixel size
Np = 20; %Number of photons per pixel

%Generating data
Data = genFLIM(EmpParam,RhoIn,D,Np);

%% Parameters used in the algorithm

PixNum = 10;
BNP.PerSample = 200;
BNP.M = M;
BNP.Alpha = 1; %Shape parameter of gamma prior on inverse of lifetimes (lambda)
BNP.Beta = 5; %Scale parameter of gamma prior on lambda
BNP.Alpha_Prop = 2000; %Parameter of proposal distribution of lambda
BNP.NJump = 200000; %Number of samples (iterations)
BNP.D = EmpParam.PixelSize/2; %Grid size in Gaussian process (pixel)
BNP.T = 1; %GP prior parameter
BNP.L = max(PixNum+2*EmpParam.OmegaX)*EmpParam.PixelSize*0.25; %GP prior parameter
BNP.N = 5; %Photons can be detected up to N pulses after excitation pulse
BNP.Sig_GP = 0.005;
BNP.Alpha_Rho = 1000;
BNP.Sig_Xi = 0.5;
BNP.Sig_Prior_Xi = 3;
BNP.Gamma = 3;
BNP.DEBUG = 0;

%% Making Inference

Data=reshape(Data,[YPixNum XPixNum]);
for ii = 1:size(Data,1)
    for jj = 1:size(Data,2)
        Data(ii,jj).X_Confocal = (0.5+(jj-1))*EmpParam.PixelSize; 
        Data(ii,jj).Y_Confocal = (0.5+(ii-1))*EmpParam.PixelSize;
    end
end
Data = Data(:);
Lambda_init = gamrnd(1,2,[1 M]);

for mm = 1:M
    Rho_init(mm).Rho = 0.5+10*rand()*ones(size(Data));
end
EmpParam.Mu = ones(1,M);
Str = 'SimData';

tic();
Chain=runBNPs_FLIM(Data,EmpParam,BNP,Rho_init,Lambda_init,Str);
T = toc();
fprintf('It took %f s to analyze this data.\n',T)

save(Str,'Chain','-v7.3')

%%

Tau = zeros(3000,6);
Load = zeros(3000,1);
for ii = 1:3000
    Load(ii) = sum(Chain(ii).Loads);
    Tau(ii,:) = 1./Chain(ii).Lambda;
end

XX = 1:3000;
figure;plot(Load,'linewidth',1.2);ylim([0 6])
hold;plot(XX,3*ones(size(XX)),'--r','linewidth',1.1)
xlabel('#samples');ylabel('active loads')

figure;plot(Tau);ylim([0 7])
hold;plot(XX,ones(size(XX)),'r--','linewidth',1.5)
plot(XX,1.8*ones(size(XX)),'r--','linewidth',1.5)
plot(XX,4.5*ones(size(XX)),'r--','linewidth',1.5)
xlabel('#samples');ylabel('lifetime (ns)')
set(gca,'Ytick',[1 1.8 4.5])
