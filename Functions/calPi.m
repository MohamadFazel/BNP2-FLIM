function [P,P0A] = calPi(EmpParam,RhoIn,Loads,TestP,D,ii,IDTestP,NormalizeFlag)
%calPi finds the probability of photons from different species given concentrations
%
%INPUT:
%   EmpParam: Structure containing parameters of the experiment
%   RhoIn: Structure containing inline functions of concentrations
%   ii: Index of the confocal spot
%
%OUTPUT:
%   P: Vector of calculated probabilities
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2022)
%   

M = length(RhoIn);
P0A = zeros(M,1);
Pi_At = zeros(M,1);
for mm = 1:M
    if Loads(mm)
        IntA = discreteInt(TestP,RhoIn(mm).Rho,EmpParam,D,ii,IDTestP);
        P0A(mm) = exp(-EmpParam.Mu(mm)*EmpParam.Dp*IntA);
    else
        P0A(mm) = 1;
    end
end
%P0A = P0A(P0A~=1);
Mm = length(P0A);
for mm = 1:Mm
    Pi_At(mm) = (1-P0A(mm))*prod(P0A)/P0A(mm);
end

PSum = sum(Pi_At);
if PSum == 0; PSum = 1; end
P=zeros(Mm,1);
if NormalizeFlag
    for mm = 1:Mm
        P(mm) = Pi_At(mm)/PSum;
    end
else
    for mm = 1:Mm
        P(mm) = Pi_At(mm);
    end
end

end
