function [Rho_IndP,Rho_Found,LikeL,GPprior]=sampleRho_Eliptical(Data,EmpParam,Chol_IndP,...
    Lambda,Rho_IndP,Rho_Found,TestP,IDTestP,Matrix,BNP,MeanXi,Loads)
%sampleRho_Eliptical uses the eliptical sampling algorithm to sample the
%lifetime maps
%
%INPUT:
%   Data: Structure array containing photon arrival times (ns)
%   EmpParam: Structure containing parameters of the experiment
%   Chol_IndP: Choleski matrix corresponding to the GP covariance matrix  
%   Lambda: inverse lifetimes
%   Rho_IndP: Lifetime maps at inducing points (pixel centers)
%   Rho_Found: Lifetime maps at the desired grid-points (test points)
%   TestP: Set of test points 
%
%OUTPUT:
%   Rho_IndP: Sampled lifeime maps at inducing points
%   Rho_Found: Sampled lifetime maps at test points
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2022) 
%

%Updating first Rho
Rho_IndPtmp = Rho_IndP;
Rho_tmp = Rho_Found;

M = BNP.M;
Ntmp = BNP.N;
D = BNP.D;
DEBUG = BNP.DEBUG;

%proposing new maps using the prior
Theta_Interval1 = [0,2*pi];
Theta1 = (Theta_Interval1(2) - Theta_Interval1(1))*rand() + Theta_Interval1(1);
Theta_Interval1(1) = Theta1 - 2*pi;
Theta_Interval1(2) = Theta1;

Theta_Interval2 = [0,2*pi];
Theta2 = (Theta_Interval2(2) - Theta_Interval2(1))*rand() + Theta_Interval2(1);
Theta_Interval2(1) = Theta2 - 2*pi;
Theta_Interval2(2) = Theta2;

u = log(rand());
%Accept = 0;
BreakWhile = 0;
NModif = 2;
mm = randperm(M,NModif);

NormFlag = 1;
LogLike_Old = 0;
if sum(Loads(mm)) ~= 0
    for ii = 1:length(Data)
        LikeT = calLikelihood(Data(ii).Dt,Lambda,Ntmp,EmpParam.T_IRF,EmpParam.T,EmpParam.Sig_IRF);
        [Pi_S,P0A] = calPi(EmpParam,Rho_Found,Loads,TestP,D,ii,IDTestP,NormFlag);
        Pi_S = repmat(Pi_S',[length(Data(ii).Dt),1]);
        LikeS = sum(Pi_S.*LikeT,2);
        LogLike_Old = LogLike_Old + sum(log(LikeS));
        LogLike_Old = LogLike_Old + ...
           sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
           sum(Data(ii).W==0)*(sum(log(P0A)));
    end
end

Xi1 = log(Rho_IndP(mm(1)).Rho);
tXi1 = MeanXi(mm(1)) + ((randn([1,length(Rho_IndP(1).Rho)])*Chol_IndP)');
Xi2 = log(Rho_IndP(mm(2)).Rho);
tXi2 = MeanXi(mm(2)) + ((randn([1,length(Rho_IndP(1).Rho)])*Chol_IndP)');

while BreakWhile == 0
    
    if Theta_Interval1(2) - Theta_Interval1(1) < 0.002 ...
            && Theta_Interval2(2) - Theta_Interval2(1) < 0.002
       BreakWhile = 1;
    end
    

    tXi_Slice1 = Xi1*cos(Theta1)+tXi1*sin(Theta1);
    Rho_IndPtmp(mm(1)).Rho = exp(tXi_Slice1);
    Rho_tmp(mm(1)).Rho = exp(Matrix*tXi_Slice1);

    tXi_Slice2 = Xi2*cos(Theta2)+tXi2*sin(Theta2);
    Rho_IndPtmp(mm(2)).Rho = exp(tXi_Slice2);
    Rho_tmp(mm(2)).Rho = exp(Matrix*tXi_Slice2);
    
    LogLike_Prop = 0;
    %log likelihood of the proposed params
    if sum(Loads(mm)) ~= 0
        for ii = 1:length(Data)
             LikeT = calLikelihood(Data(ii).Dt,Lambda,Ntmp,EmpParam.T_IRF,EmpParam.T,EmpParam.Sig_IRF);
             [Pi_S,P0A] = calPi(EmpParam,Rho_tmp,Loads,TestP,D,ii,IDTestP,NormFlag);
             Pi_S = repmat(Pi_S',[length(Data(ii).Dt),1]);
             LikeS = sum(Pi_S.*LikeT,2);
             LogLike_Prop = LogLike_Prop + sum(log(LikeS));
             LogLike_Prop = LogLike_Prop + ...
                 sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
                 sum(Data(ii).W==0)*(sum(log(P0A)));
        end
    end
    
    if (LogLike_Prop - LogLike_Old) > u
         Rho_Found = Rho_tmp; 
         Rho_IndP = Rho_IndPtmp;
         BreakWhile = 1;
    end
    if Theta1 < 0
        Theta_Interval1(1) = Theta1;
    else
        Theta_Interval1(2) = Theta1; 
    end
    Theta1 = (Theta_Interval1(2) - Theta_Interval1(1))*rand() + Theta_Interval1(1);
    
    if Theta2 < 0
        Theta_Interval2(1) = Theta2;
    else
        Theta_Interval2(2) = Theta2; 
    end
    Theta2 = (Theta_Interval2(2) - Theta_Interval2(1))*rand() + Theta_Interval2(1);
    
end

%log prior ratio for b=1
LogPriorProp = 0;
LogPriorOld = 0;
for mm = 1:M
    LogPriorProp = LogPriorProp -0.5*sum(((log(Rho_IndPtmp(mm).Rho)-(MeanXi(mm)))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndPtmp(mm).Rho)*log(2*pi)/2;
    LogPriorOld = LogPriorOld -0.5*sum(((log(Rho_IndP(mm).Rho)-(MeanXi(mm)))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndP(mm).Rho)*log(2*pi)/2;
end

LikeL = LogLike_Prop;
GPprior = LogPriorProp;

for mm = 1:M
    Rho_Found(mm).Rho = single(Rho_Found(mm).Rho);
end

if DEBUG
    figure(111);
    L = [max(EmpParam.Xi_Y),max(EmpParam.Xi_X)]/EmpParam.PixelSize + 0.5;
    [Xg,Yg] = meshgrid((0.5:L(2)-0.5),(0.5:L(1)-0.5));
    for mm = 1:M
        Rho = Rho_IndP(mm).Rho;
        Rho1 = reshape(Rho,L);
        if Loads(mm) == 1
            surf(Xg,Yg,Rho1)
        end
        if mm == 1
            hold(gca,'on');
        end
    end
    xlabel('X');ylabel('Y')
    pause(0.02)
    hold off;
end

end
