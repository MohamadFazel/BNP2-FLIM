function findIRF(Dt,Lifetime,T,N,Bg,NChain,Alpha)
%findIRF find IRF parameters using calibration data with a single species.
%The IRF assumes to be Gaussian and a single lifetime component with
%exponential lifetime decay. The background is assumed to be a uniform
%distribution with not IRF effect.
%
%INPUT:
%   Dt: Photon arrival times (ns)
%   Lifetime: Mean lifetime of the species in the calibration data
%   T: pulse period (time preiod between 2 consecutive pulses) (ns)
%   N: Number of pretiouse pulse to be considered
%   Bg: Binary parameter indicating if there is background or not (Defalt: 0)
%   NChain: Number of samples to be taken (Default: 10000)
%   Alpha: Parameter of the gamma proposal distribution (Default: 1000)
%
%OUTPUT:
%   None
%
%Created by:
%   Mohamadreza Fazel (Presse lab, 2021)
%
if nargin < 5
    Bg = 0; 
end
if nargin < 6
    NChain = 10000;
end
if nargin < 7
    Alpha = 5000;
end

%Chains
Tau_Chain = zeros(NChain,1);
Sig_Chain = zeros(NChain,1);
W_Chain = zeros(NChain,1);

%Initializing the chain
Tau_Chain(1) = 12;
Sig_Chain(1) = 1;
if Bg
    W_Chain(1) = 0.1;
else
    W_Chain(1) = 0;
end

%inline functions
Lambda = 1/Lifetime;
LExp = @(Delt,tau,sig,n) (Lambda/2)*exp((Lambda/2)*(2*(tau-Delt-n*T)+Lambda*sig^2));
LErf = @(Delt,tau,sig,n) erfc((tau-Delt-n*T+Lambda*sig^2)/(sig*sqrt(2)));
ExtraTerm = @(Delt,tau,Sig) 1/T;

for ii = 2:NChain
    
    Tau_Current = Tau_Chain(ii-1);
    Sig_Current = Sig_Chain(ii-1);
    W_Current = W_Chain(ii-1);
    
    %Sampling Tau
    Tau_Prop = gamrnd(Alpha,Tau_Current/Alpha);
  
    Ltmp_Current1 = 0;
    Ltmp_Prop1 = 0;
    for nn = 0:N
        Ltmp_Current1 = Ltmp_Current1 + LExp(Dt,Tau_Current,Sig_Current,nn)...
            .*LErf(Dt,Tau_Current,Sig_Current,nn); 
        Ltmp_Prop1 = Ltmp_Prop1 + LExp(Dt,Tau_Prop,Sig_Current,nn)...
         .*LErf(Dt,Tau_Prop,Sig_Current,nn);   
    end
    
    Ltmp_Current2 = (1-W_Current)*Ltmp_Current1 + W_Current*ExtraTerm(Dt,Tau_Current,Sig_Current);
    Ltmp_Prop2 = (1-W_Current)*Ltmp_Prop1 + W_Current*ExtraTerm(Dt,Tau_Prop,Sig_Current);
    
    DLogL = sum(log(Ltmp_Prop2) - log(Ltmp_Current2));
    DPrior = sum(log(gampdf(Tau_Prop,3,4)) - log(gampdf(Tau_Current,3,4)));
    DProp = sum(log(gampdf(Tau_Current,Alpha,Tau_Prop/Alpha)) ...
        - log(gampdf(Tau_Prop,Alpha,Tau_Current/Alpha)));
    
    if DLogL + DPrior + DProp > log(rand())
        Tau_Current = Tau_Prop;
    end
    
    %Sampling Sigma
    Sig_Prop = gamrnd(Alpha,Sig_Current/Alpha);
    
    Ltmp_Current1 = 0;
    Ltmp_Prop1 = 0;
    for nn = 0:N
        Ltmp_Current1 = Ltmp_Current1 + LExp(Dt,Tau_Current,Sig_Current,nn)...
            .*LErf(Dt,Tau_Current,Sig_Current,nn); 
        Ltmp_Prop1 = Ltmp_Prop1 + LExp(Dt,Tau_Current,Sig_Prop,nn)...
         .*LErf(Dt,Tau_Current,Sig_Prop,nn);   
    end
    
    Ltmp_Current2 = (1-W_Current)*Ltmp_Current1 + W_Current*ExtraTerm(Dt,Tau_Current,Sig_Current);
    Ltmp_Prop2 = (1-W_Current)*Ltmp_Prop1 + W_Current*ExtraTerm(Dt,Tau_Current,Sig_Prop);
    
    DLogL = sum(log(Ltmp_Prop2) - log(Ltmp_Current2));
    DPrior = sum(log(gampdf(Sig_Prop,1,4)) - log(gampdf(Sig_Current,1,4)));
    DProp = sum(log(gampdf(Sig_Current,Alpha,Sig_Prop/Alpha)) ...
        - log(gampdf(Sig_Prop,Alpha,Sig_Current/Alpha)));
    
    if DLogL + DPrior + DProp > log(rand())
        Sig_Current = Sig_Prop;
    end
    
    %Sampling Weight
    if Bg 
        
        W_Prop = gamrnd(Alpha*10,W_Current/Alpha/10);

        Ltmp_Current1 = 0;
        for nn = 0:N
            Ltmp_Current1 = Ltmp_Current1 + LExp(Dt,Tau_Current,Sig_Current,nn)...
                .*LErf(Dt,Tau_Current,Sig_Current,nn);    
        end

        Ltmp_Current2 = (1-W_Current)*Ltmp_Current1 + W_Current*ExtraTerm(Dt,Tau_Current,Sig_Current);
        Ltmp_Prop2 = (1-W_Prop)*Ltmp_Current1 + W_Prop*ExtraTerm(Dt,Tau_Current,Sig_Current);

        DLogL = sum(log(Ltmp_Prop2) - log(Ltmp_Current2));
        DPrior = sum(log(betapdf(W_Prop,1,1)) - log(betapdf(W_Current,1,4)));
        DProp = sum(log(gampdf(W_Current,Alpha*10,W_Prop/Alpha/10)) ...
            - log(gampdf(W_Prop,Alpha*10,W_Current/Alpha/10)));

        if DLogL + DPrior + DProp > log(rand())
            W_Current = W_Prop;
        end
    
    end
    
    Tau_Chain(ii) = Tau_Current;
    Sig_Chain(ii) = Sig_Current;
    W_Chain(ii) = W_Current;

end

figure;plot(Tau_Chain);ylabel('IRF offset (ns)')
Tau_mean=mean(Tau_Chain(NChain-1000:NChain));

figure;plot(Sig_Chain);ylabel('IRF width (ns)')
Sig_mean=mean(Sig_Chain(NChain-1000:NChain));

figure;plot(W_Chain);ylabel('\pi')
A_mean=mean(W_Chain(NChain-1000:NChain));

DelT = 0:0.02:T+0.1;
figure;histogram(Dt,256,'normalization','pdf')
Y=0;
for nn = 0:N
    Y = Y + LExp(DelT,Tau_mean,Sig_mean,nn).*LErf(DelT,Tau_mean,Sig_mean,nn);% + ...
            %  A_mean*Extra(DelT,Tau_mean,Sig_mean,nn);
end
Y = (1-A_mean)*Y + A_mean/T;
hold;plot(DelT,Y,'r','linewidth',1.8)
xlabel('t(ns)','FontSize',16);
ylabel('PDF','FontSize',16)
text(4,0.15,'Calibration data')
text(4,0.135,sprintf('IRF offset: %f ns',Tau_mean));
text(4,0.12,sprintf('IRF std: %f ns',Sig_mean));
text(4,0.105,sprintf('Bg fraction: %f ns',A_mean));
legend('Data','Fit','location','North')

end
