function Chain=runBNPs_FLIM(Data,EmpParam,BNP,Rho_init,Lambda_init,Str,initLoads,initXi)
%runBNPs_FLIM implements complete BNPs analysis of FLIM data.
%
%INPUT:
%   Data: Structure array containing photons arrival times (ns)
%   EmpParams: Structure containing experimental parameters
%   RhoIn:
%   NJump: Number of MCMC jumps
%   Alpha: Shape parameter of gamma prior on lambda
%   Beta:  Rate parameter of gamma prior on lambda
%   Alpha_Prop: Shape parameter for lambda proposal distribution
%   S_init: Initial values fot tag parameter S (optional)
%   Lambda_init: Initial values for labmda (optional)
%
%OUTPUT:
%   Chain: Structure array containing lambda and S parameters of MCMC jumps
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2022)
%

NJump = BNP.NJump;
Alpha = BNP.Alpha;
Beta = BNP.Beta;
Alpha_Prop = BNP.Alpha_Prop;
D = BNP.D;
T = BNP.T;
L = BNP.L;
M = BNP.M;
if isfield(BNP,'Exten')
    Exten = BNP.Exten;
else
    Exten = 1;
end
if isfield(BNP,'PerSample')
    PerSample = BNP.PerSample;
else
    PerSample = 50;
end

%Initializing the chain
Chain(floor(NJump/PerSample)).Lambda = zeros(1,BNP.M);
Chain(floor(NJump/PerSample)).Loads = zeros(1,BNP.M);
Chain(floor(NJump/PerSample)).Xi = zeros([1,M]);
Chain(floor(NJump/PerSample)).LogLikelihood = 0;
Chain(floor(NJump/PerSample)).LogPosterior = 0;
Chain(floor(NJump/PerSample)).LogRhoPrior = 0;
Chain(floor(NJump/PerSample)).LogXiPrior = 0;
Chain(floor(NJump/PerSample)).LogLamPrior = 0;
Chain(floor(NJump/PerSample)).LogLoadPrior = 0;
if nargin < 7
    Chain(1).Loads = round(rand(1,BNP.M));
    Chain(1).Loads(randi(BNP.M)) = 1; 
else
    Chain(1).Loads = initLoads;
end
if nargin < 8
    Chain(1).Xi = 0.5+10*rand([1,M]);
else
    Chain(1).Xi = initXi;
end

Data = Data(:);
Cent_Confocal = zeros(length(Data),3);
for ii = 1:length(Data)
    Cent_Confocal(ii,1) = Data(ii).X_Confocal;
    Cent_Confocal(ii,2) = Data(ii).Y_Confocal;
    Cent_Confocal(ii,3) = Data(ii).Z_Confocal;
end

%Initializing the chain
if nargin >= 5
    if size(Lambda_init,1)>1
        Chain(1).Lambda = Lambda_init';
    else
        Chain(1).Lambda = Lambda_init;
    end
    Chain(1).Rho_IndP = Rho_init;
elseif nargin == 4
    Chain(1).Lambda = rand([1 M]);
    Chain(1).Rho_IndP = Rho_init;
elseif nargin == 3
    Chain(1).Lambda = rand([1 M]);
    for mm = 1:M
        Chain(1).Rho_IndP(mm).Rho = 20*rand()*ones([length(Cent_Confocal),1]);
    end
end

%Covariance matrix between inducing points
K_IndP = calCOV(Cent_Confocal,Cent_Confocal,T,L);
MinRange = -(Exten*EmpParam.PixelSize-D/2)*[1, 1]; 
MaxRange(1) = (Data(end).X_Confocal+Exten*EmpParam.PixelSize);
MaxRange(2) = (Data(end).Y_Confocal+Exten*EmpParam.PixelSize);
[Xg,Yg,Zg] = meshgrid(MinRange(1):D:MaxRange(1),MinRange(2):D:MaxRange(2),0);

TestP = [Xg(:),Yg(:),Zg(:)];
IDTestP = zeros([size(TestP,1),length(Data)]);
for ii = 1:size(Cent_Confocal,1)
     Dis = pdist2(TestP(:,1:2),Cent_Confocal(ii,1:2));
     IDTestP(Dis<3*EmpParam.OmegaX,ii) = ii;
end

%Covariance matrix between inducing and test points
K_test_IndP = calCOV(TestP,Cent_Confocal,T,L);
%Choleski decomposition of covariance matrix
Chol_IndP = cholcov(K_IndP+1000*eps*eye(size(K_IndP)));
%The intermediate matrix that relate indusing points to test points
Matrix = K_test_IndP/K_IndP;
%Initializing the profiles at the test points based on the given values
%for test points and the correlation matrix.
for mm = 1:M
    Chain(1).Rho(mm).Rho = exp(Matrix*log(Chain(1).Rho_IndP(mm).Rho)); 
end

Ind = 1;
ThisRho_IndP = Chain(1).Rho_IndP;
ThisRho = Chain(1).Rho;
ThisLambda = Chain(1).Lambda;
ThisXi = Chain(1).Xi;
ThisLoads = Chain(1).Loads;
AcceptXi = 0;
AcceptLam = 0;
for jj = 2:NJump
    
    if jj/10000 == floor(jj/10000)
        fprintf('Jump %d out of %d\n',jj,NJump);
        fprintf('Accepted Lifetime jumps: %d\n',AcceptLam);
        fprintf('Accepted Xi jumps: %d\n',AcceptXi);
    end
    
    [ThisRho_IndP,ThisRho,~,GPprior]= ...
        sampleRho_Eliptical(Data,EmpParam,Chol_IndP,ThisLambda,ThisRho_IndP,...
        ThisRho,TestP,IDTestP,Matrix,BNP,ThisXi,ThisLoads);
    
    %Sampling means of GPs
    [ThisXi,AcceptXi,XiPrior] = sampleXi(ThisRho_IndP,ThisXi,BNP,Chol_IndP,ThisLoads,AcceptXi);
         
    %Sampling lifetimes
    [ThisLambda,~,AcceptLam,LamPrior] = sampleLambda(Data,EmpParam,ThisRho,...
        ThisLambda,Alpha,Beta,Alpha_Prop,BNP.N,TestP,IDTestP,D,ThisLoads,AcceptLam);
    
    %Sampling loads via beta-Bernoulli process
    [ThisLoads,Like,LoadPrior] = sampleSingleLoads(Data,ThisLambda,ThisLoads,ThisRho,...
         EmpParam,BNP,TestP,IDTestP); 

    %Save every "PerSample" step of the chain 
    if jj/PerSample == floor(jj/PerSample)
        Ind = Ind + 1;
        Chain(Ind).Rho_IndP = ThisRho_IndP;
        Chain(Ind).Rho = ThisRho;
        Chain(Ind).Lambda = ThisLambda;
        Chain(Ind).Xi = ThisXi;
        Chain(Ind).Loads = ThisLoads;
        Chain(Ind).LogLikelihood = Like;
        Chain(Ind).LogRhoPrior = GPprior;
        Chain(Ind).LogXiPrior = XiPrior;
        Chain(Ind).LogLamPrior = LamPrior;
        Chain(Ind).LogLoadPrior = LoadPrior;
        Chain(Ind).LogPosterior = Like + GPprior + XiPrior + LamPrior + LoadPrior;
    end
             
    %Saving chunks of the chain
    if jj/100000 == floor(jj/100000)
        Num = (jj/100000-1);
        SInd = Num*floor(50000/PerSample)+1;
        tChain = Chain(SInd:Ind);
        save(sprintf('%s_%d',Str,Num),'tChain','-v7.3');
        sprintf('Part %d of the chain is saved\n',Num);
    end
      
end

end