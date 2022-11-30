function [Xi,AcceptXi,LogPrior] = sampleXi(Rho_IndP,Xi_Old,BNP,Chol_IndP,Loads,AcceptXi)

Sig_Xi = BNP.Sig_Xi;
Sig_Prior_Xi = BNP.Sig_Prior_Xi;
M = BNP.M;

%proposing new Xi
Xi_Prop = Xi_Old+Sig_Xi*randn([1,M]);

%log likelihood ratio
LogLikeProp = 0;
LogLikeOld = 0;
for mm = 1:M
    LogLikeProp = LogLikeProp -0.5*sum(((log(Rho_IndP(mm).Rho)-((Xi_Prop(mm))))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndP(mm).Rho)*log(2*pi)/2;
    LogLikeOld = LogLikeOld -0.5*sum(((log(Rho_IndP(mm).Rho)-((Xi_Old(mm))))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndP(mm).Rho)*log(2*pi)/2;
end
LogLikeRatio = LogLikeProp - LogLikeOld;

%log prior ratio
LogPriorProp = log(normpdf(Xi_Prop,-1,Sig_Prior_Xi));
LogPriorOld = log(normpdf(Xi_Old,-1,Sig_Prior_Xi));
LogPriorRatio = sum(LogPriorProp-LogPriorOld);

if LogLikeRatio + LogPriorRatio > log(rand())
    Xi = Xi_Prop;
    Xi(Loads==0) = Sig_Prior_Xi*randn(1,sum(Loads==0));
    AcceptXi = AcceptXi + 1;
    LogPrior = sum(LogPriorProp);
else
    Xi = Xi_Old;
    LogPrior = sum(LogPriorOld);
end

end
