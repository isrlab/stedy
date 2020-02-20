function Ksys = genKsys(x,tData)
% Function to generate Ksys from section 3.2. At equilibrium.

ns = numel(tData.N);
q = x(1:ns);
Ksys = zeros(ns,ns);

Y = tData.Y;
X = tData.X;
for k = 1:tData.nStr
    bk = tData.listY{k}*q;
    lsk = norm(bk);
    Kstr1 = tData.sigmaEq(k)*Y{k};
    Kstr2 = tData.K(k)*tData.Lk(k)/(lsk^3)*Y{k}*(q*q.')*Y{k};
    Ksys = Ksys + Kstr1 + Kstr2;
end
if(tData.Compressible)
    for j = 1:tData.nBar
        bj = tData.listX{j}*q;
        lbj = norm(bj);
        Kbar1 = tData.psiEq(j)*X{j};
        Kbar2 = tData.bars.listK(j)*tData.bars.L0(j)/(lbj^3)*X{j}*(q*q.')*X{j};
        Ksys = Ksys + Kbar1 + Kbar2;
    end
end