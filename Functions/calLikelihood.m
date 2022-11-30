function LikeOut = calLikelihood(Data,Lambda,Ntmp,T_IRF,T,Sig_IRF)
%calLikelihood finds the likelihood of the proposed lamda
%
%INPUT:
%   Data: Structure array containing photon arrival times (ns)
%   EmpParam: Structure containing parameters of the experiment
%   S:  Sampled S in the previous iteration
%   Lambda: Proposed lambda (1/ns)
%
%OUTPUT:
%   Likelihood: Calculated likelihood of the proposed lambda
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2022)
%

if length(Lambda)==1
    LambdaS = Lambda*ones(length(Data),1);
else
    LambdaS = repmat(Lambda,[length(Data),1]);
end

LambdaT = repmat(LambdaS,[1,1,Ntmp+1]);
DataT = repmat(Data,[1,length(Lambda),Ntmp+1]);
tNT = zeros(1,1,Ntmp+1);
tNT(1,1,:) = (0:Ntmp);
NT = repmat(tNT,[length(Data),length(Lambda),1]);

LikeExp = (LambdaT/2).*exp((LambdaT/2).*(2*(T_IRF-DataT-NT*T) + ...
        LambdaT*Sig_IRF.^2));
LikeErf = erfc((T_IRF-DataT-NT*T+LambdaT*Sig_IRF.^2) ...
        /(sqrt(2)*Sig_IRF)); 
    
LikeOut = sum(LikeExp.*LikeErf,3);
LikeOut(isnan(LikeOut)) = 0;
LikeOut(isinf(LikeOut)) = 0;
      
end
