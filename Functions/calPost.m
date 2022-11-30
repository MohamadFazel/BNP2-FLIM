function [LogP,Like_Dt] = calPost(Data,Chain,S,EmpParam,BNP,Ntmp,TestP,IDTestP,Chol_IndP)

Like_Dt = 0;
D = BNP.D;
Lambda = Chain.Lambda;
%log Likelihood Dt
for ii = 1:length(Data)
     Stmp = S(ii).S(Data(ii).W==1);
     Like_Dt = Like_Dt + calLikelihood(Data(ii).Dt,EmpParam,Stmp,Lambda,Ntmp,...
         EmpParam.T_IRF,EmpParam.T,EmpParam.Sig_IRF);
end

%LogLikelihood W,S
NormFlag = 1;
Like_Rho = 0;
for ii = 1:length(Data)
    [Pi_S,P0A] = calPi(EmpParam,Chain.Rho,TestP,D,ii,IDTestP,NormFlag);
     St = S(ii).S;
     Like_Rho = Like_Rho + sum(log(Pi_S(St))) + ...
        sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
        sum(Data(ii).W==0)*(sum(log(P0A)));
end

%Prior Lambda
PriorGamma = sum(log(gampdf(Lambda,BNP.Alpha,BNP.Beta)));
%Prior Rho
PriorRho = 0;
for ii = 1:length(Chain.Rho)
    tRho = Chain.Rho_IndP(ii).Rho';
    d = length(tRho);
    logSqrtDetSigma = sum(log(diag(Chol_IndP)));
    xRinv = (log(tRho)-Chain.Xi(ii)) / Chol_IndP;
    quadform = sum(xRinv.^2,2);
%     [~,logY]=mvnpdf_Log(log(tRho),Chain.Xi(ii),Chol_IndP*Chol_IndP');
    PriorRho = PriorRho + sum(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);
end

%Prior \nu
PriorXi = sum(log(normpdf(Chain.Xi,0,50)));

LogP = Like_Dt + Like_Rho + PriorGamma + PriorRho + PriorXi;

end

