function [Probs,LikeOut] = calProbs(Data,EmpParam,Lambda,Pi_T,Ntmp)
%calProbs finds probabilities used in categorical distribution to sample S
%
%INPUT:
%   Data: Structure array containing photon arrival times (ns)
%   EmpParam: Structure containing parameters of the experiment
%   Lambda: Sampled lambda in the previous iteration (1/ns)
%   Pi_T: Relative probability of photons coming from different species
%
%OUTPUT:
%   Probs: The probabilities of a given photon coming from different species
%
%Created by:
% Mohamadreza Fazel (Presse Lab, 2022)
%

LambdaS = repmat(Lambda,[length(Data),1,Ntmp(end)+1]);
DataS = repmat(Data,[1,length(Lambda),Ntmp(end)+1]);
Nt = 0:Ntmp(end);
Nt = reshape(Nt,[1,1,length(Nt)]);
NS = repmat(Nt,[length(Data),length(Lambda),1]);
LikeExp = (LambdaS/2).*exp((LambdaS/2).*(2*(EmpParam.T_IRF-(DataS+NS*EmpParam.T)) + ...
        LambdaS*(EmpParam.Sig_IRF^2))); 
LikeErf = erfc((EmpParam.T_IRF-(DataS+NS*EmpParam.T)+LambdaS*EmpParam.Sig_IRF.^2) ...
        /(sqrt(2)*EmpParam.Sig_IRF));    
LikeOut = sum(LikeExp.*LikeErf,3); 
        
PiT = repmat(Pi_T,[length(Data),1]);
Probs1 = LikeOut.*PiT;
NScale = repmat(sum(Probs1,2),[1,length(Lambda)]); 
Probs = Probs1./NScale;

end