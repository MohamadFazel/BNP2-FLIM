%Example script for simulated FLIM data similar to those in Fig. 2 in the 
%main text and Supplementary Fig. 1 with lifetimes of 1ns, 1.8ns and 4.5ns.
%It will take ~2 hours to run this cript.
addpath('Functions')
%Load data
load('Data_20Ph_per_pixel.mat');

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

M=6; %Number of loads used in our beta-Bernoulli process 
EmpParam.Xi_X = [];
XPixNum = 20; %Number of pixels along the X-axis
YPixNum = 5; %Number of pixels along the Y-axis
EmpParam.PixelSize = 0.3922; %data pixel size (mu)
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

%% Parameters used in the algorithm

BNP.PerSample = 100; %save every 100 samples within the output chain
BNP.M = M; %Number of loads used in the beta-Bernoulli process
BNP.Alpha = 1; %Shape parameter of gamma prior on inverse of lifetimes (lambda)
BNP.Beta = 5; %Scale parameter of gamma prior on lambda
BNP.Alpha_Prop = 2000; %Parameter of proposal distribution of lambda
BNP.NJump = 50000; %Number of samples (iterations)
BNP.D = EmpParam.PixelSize/2; %Grid size in Gaussian process (pixel)
BNP.T = 1; %GP prior parameter
BNP.L = 1; %GP prior parameter (mu)
BNP.N = 5; %Photons can be detected up to N pulses after excitation pulse
BNP.Sig_GP = 0.005; %parameter of proposal distribution for GP mean
BNP.Alpha_Rho = 1000; %parameter of prior on Rho
BNP.Sig_Xi = 0.5; %parameter of proposal distribution
BNP.Sig_Prior_Xi = 3; %parameter of prior
BNP.Gamma = 3; %Expected number of species 
BNP.DEBUG = 0;

%% Making Inference

%initialize lambda
Lambda_init = gamrnd(1,2,[1 M]);

%initialize lifetime maps
for mm = 1:M
    Rho_init(mm).Rho = 0.5+10*rand()*ones(size(Data));
end
EmpParam.Mu = ones(1,M);
%The output is saved in a file with the following name
Str = 'SimData';

tic();
Chain=runBNPs_FLIM(Data,EmpParam,BNP,Rho_init,Lambda_init,Str);
T = toc();
fprintf('It took %f s to analyze this data.\n',T)
save(Str,'Chain','-v7.3')

%% Displaying results: plot of number of found species and histogram of lifetimes

Tau = zeros(length(Chain),BNP.M);
Load = zeros(length(Chain),BNP.M);
for ii = 1:length(Chain)
    Load(ii,:) = Chain(ii).Loads;
    Tau(ii,:) = 1./Chain(ii).Lambda;
end

XX = 1:length(Chain);
figure;plot(sum(Load,2),'linewidth',1.2);ylim([0 6])
hold;plot(XX,3*ones(size(XX)),'--r','linewidth',1.1)
xlabel('#samples');ylabel('active loads');xlim([0 length(Chain)])

tTau = Tau(100:end,:);
tLoad = Load(100:end,:);
figure;HistHandle=histogram(tTau(tLoad==1),'BinWidth',0.25,'normalization','pdf');
Xx = 0:0.1:HistHandle.Parent.YLim(2);
hold;plot(ones(size(Xx)),Xx,'r--','linewidth',1.5)
plot(1.8*ones(size(Xx)),Xx,'r--','linewidth',1.5)
plot(4.5*ones(size(Xx)),Xx,'r--','linewidth',1.5)
ylabel('pdf');xlabel('lifetime (ns)')
set(gca,'Xtick',[1 1.8 4.5])
