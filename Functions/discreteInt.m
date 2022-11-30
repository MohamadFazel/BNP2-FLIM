function Integ = discreteInt(TestP,Rho,EmpParam,D,ii,IDTestP)
Xp = EmpParam.Xi_X(ii);
Yp = EmpParam.Xi_Y(ii);
Zp = EmpParam.Xi_Z(ii);
if nargin < 6 || isempty(IDTestP) 
    Dis = pdist2(TestP,[Xp,Yp,Zp]);
    ID = Dis<3*EmpParam.OmegaX;
else
   ID = IDTestP(:,ii)~=0;  
end
TestP_tmp = TestP(ID,:);

Integ=sum(Rho(ID).*EmpParam.PSF(TestP_tmp(:,1),TestP_tmp(:,2),Xp,Yp,...
    EmpParam.OmegaX,EmpParam.OmegaZ))*D^2;

end

