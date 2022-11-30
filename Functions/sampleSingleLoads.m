function [Loads,LikeL,LoadPrior] = sampleSingleLoads(Data,Lambda,Loads,...
             Rho_Found,EmpParam,BNP,TestP,IDTestP)

D=BNP.D;
M = BNP.M;
Ntmp = BNP.N;

ID = randperm(M,2);
tLoads = Loads;
Loads1 = [0 0; 1 0; 0 1; 1 1];
PostCombin = zeros(4,1);
LikeLVec = zeros(4,1);

NormFlag = 1;
for mm = 1:4
    tLoads(ID) = Loads1(mm,:);
    for ii = 1:length(Data)
        LikeT = calLikelihood(Data(ii).Dt,Lambda,Ntmp,EmpParam.T_IRF,EmpParam.T,EmpParam.Sig_IRF);
        [Pi_S,P0A] = calPi(EmpParam,Rho_Found,tLoads,TestP,D,ii,IDTestP,NormFlag);
        Pi_S = repmat(Pi_S',[length(Data(ii).Dt),1]);
        LikeS = sum(Pi_S.*LikeT,2);
        LikeLVec(mm) = LikeLVec(mm) + sum(log(LikeS));
        PostCombin(mm) = PostCombin(mm) + sum(log(LikeS));
        PostCombin(mm) = PostCombin(mm) + ...
        sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
        sum(Data(ii).W==0)*(sum(log(P0A)));
    end
    PostCombin(mm) = PostCombin(mm) + ...
        sum(tLoads==1)*log(1/(1+(M-1)/BNP.Gamma)) ...
        + sum(tLoads==0)*log(1-1/(1+(M-1)/BNP.Gamma));
end
Ind = categSample( PostCombin);

Loads(ID) = Loads1(Ind,:);
LikeL = LikeLVec(Ind);

LoadPrior = sum(Loads==1)*log(1/(1+(M-1)/BNP.Gamma)) ...
        + sum(Loads==0)*log(1-1/(1+(M-1)/BNP.Gamma));
end

function Ind = categSample(LogP)
tInd = isinf(LogP);
if sum(LogP(tInd) > 0) ~= 0
    ptInd = tInd(LogP(tInd)>0);
    Ind = ptInd(randi(numel(ptInd)));
    return;
end
LogP = LogP-max(LogP);
P = exp(LogP)./sum(exp(LogP(~isinf(LogP)&~isnan(LogP))));
P(isnan(P)) = 0;
Ind = find(cumsum(P)-rand()>=0,1);
end
