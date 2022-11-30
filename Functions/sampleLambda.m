function [Lambda,LikeL,AcceptLam,LamPrior] = sampleLambda(Data,EmpParam,Rho,Lambda_Old,Alpha,...
    Beta,Alpha_Prop,Ntmp,TestP,IDTestP,D,Loads,AcceptLam)
%sampleLambda uses a MH algorithm to sample lambda
%
%INPUT:
%   Data: Structure array containing photon arrival times (ns)
%   EmpParam: Structure containing parameters of the experiment
%   S: Sampled S parameters in the previous iteration
%   Lambda_Old: Sampled lambda in the previous iteration
%   Alpha: Shape parameter of gamma prior on lambda
%   Beta:  Rate parameter of gamma prior on lambda
%   Alpha_Prop: Shape parameter for lambda proposal distribution
%
%OUTPUT:
%   Lambda: Sampled lambda
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2022) 
%

%loglikelihood for the current params
NormFlag = 1;
Like_Old = 0;
for ii = 1:length(Data)
     LikeT = calLikelihood(Data(ii).Dt,Lambda_Old,Ntmp,EmpParam.T_IRF,EmpParam.T,EmpParam.Sig_IRF);
     [Pi_S,~] = calPi(EmpParam,Rho,Loads,TestP,D,ii,IDTestP,NormFlag);
     Pi_S = repmat(Pi_S',[length(Data(ii).Dt),1]);
     LikeS = sum(Pi_S.*LikeT,2);
     Like_Old = Like_Old + sum(log(LikeS));
end

%proposing new values
Lambda_Prop = gamrnd(Alpha_Prop,Lambda_Old/Alpha_Prop);

%loglikelihood for the proposed values
Like_Prop = 0;
for ii = 1:length(Data)
    LikeT = calLikelihood(Data(ii).Dt,Lambda_Prop,Ntmp,EmpParam.T_IRF,EmpParam.T,EmpParam.Sig_IRF);
    [Pi_S,~] = calPi(EmpParam,Rho,Loads,TestP,D,ii,IDTestP,NormFlag);
    Pi_S = repmat(Pi_S',[length(Data(ii).Dt),1]);
    LikeS = sum(Pi_S.*LikeT,2);
    Like_Prop = Like_Prop + sum(log(LikeS));
end
LikeRatio = sum(Like_Prop-Like_Old);

%log prior ratio
LamPropPrior = log(gampdf(Lambda_Prop,Alpha,Beta));
LamOldPrior = log(gampdf(Lambda_Old,Alpha,Beta));
PriorRatio = sum(LamPropPrior - LamOldPrior);

%Proposal ratio
PropRatio = sum(log(gampdf(Lambda_Old,Alpha_Prop,Lambda_Prop/Alpha_Prop))...
     -log(gampdf(Lambda_Prop,Alpha_Prop,Lambda_Old/Alpha_Prop)));
 
A = LikeRatio+PriorRatio+PropRatio;
if A > log(rand())
    Lambda = Lambda_Prop;
    LikeL = Like_Prop;
    AcceptLam = AcceptLam + 1;
    LamPrior = sum(LamPropPrior);
else
    Lambda = Lambda_Old;
    LikeL = Like_Old;
    LamPrior = sum(LamOldPrior);
end

end
