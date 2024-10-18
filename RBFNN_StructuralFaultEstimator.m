function [wHatX,wHatY,wHatZ,wHatPhi,wHatTheta,wHatPsi, FaultHat] = ...
    RBFNN_StructuralFaultEstimator(t,e,eDot,S,bij,cij,GAMMA,...
    nKernel,WHAT_OLDx,WHAT_OLDy,WHAT_OLDz,...
    WHAT_OLDphi,WHAT_OLDtheta,WHAT_OLDpsi, bPrime,LandaPrime)
    
    n = 6;
    Ts = t(2) - t(1);
    FaultHat = zeros(n,1);
    wHatX = 0.01*ones(nKernel,1);
    wHatY = 0.01*ones(nKernel,1);
    wHatZ = 0.01*ones(nKernel,1);
    wHatPhi = 0.01*ones(nKernel,1);
    wHatTheta = 0.01*ones(nKernel,1);
    wHatPsi = 0.01*ones(nKernel,1);

    hx = wHatX;
    hy = wHatX;
    hz = wHatX;
    hPhi = wHatX;
    hTheta = wHatX;
    hPsi = wHatX;
    
    gammaX = GAMMA(1);
    gammaY = GAMMA(2);
    gammaZ = GAMMA(3);
    gammaPhi = GAMMA(4);
    gammaTheta = GAMMA(5);
    gammaPsi = GAMMA(6);
    
    xiX = [e(1)                 eDot(1)]';
    xiY = [e(2)                 eDot(2)]';
    xiZ = [e(3)                 eDot(3)]';
    xiPhi = [e(4)        	    eDot(4)]';
    xiTheta = [e(5)           eDot(5)]';
    xiPsi = [e(6)               eDot(6)]';
    
    sx = S(1);
    sy = S(2);
    sz = S(3);
    sPhi = S(4);
    sTheta = S(5);
    sPsi = S(6);
    
    for j=1:nKernel
        
        hx(j) = exp(-norm(xiX-cij(1,j))^2/(2*(bij(1,j))^2));
        wHatX(j) = WHAT_OLDx(j)+...
        Ts*(1/gammaX)*bPrime(1)*LandaPrime*abs(eDot(1))^(LandaPrime-1)*sx*hx(j);   % Solve wHatDot
        
        hy(j) = exp(-norm(xiY-cij(2,j))^2/(2*(bij(2,j))^2));
        wHatY(j) = WHAT_OLDy(j)+...
        Ts*(1/gammaY)*bPrime(2)*LandaPrime*abs(eDot(2))^(LandaPrime-1)*sy*hy(j);
        
        hz(j) = exp(-norm(xiZ-cij(3,j))^2/(2*(bij(3,j))^2));
        wHatZ(j) = WHAT_OLDz(j)+...
        Ts*(1/gammaZ)*bPrime(3)*LandaPrime*abs(eDot(3))^(LandaPrime-1)*sz*hz(j);
        
        hPhi(j) = exp(-norm(xiPhi-cij(4,j))^2/(2*(bij(4,j))^2));
        wHatPhi(j) = WHAT_OLDphi(j)+...
        Ts*(1/gammaPhi)*bPrime(4)*LandaPrime*abs(eDot(4))^(LandaPrime-1)*sPhi*hPhi(j);

        hTheta(j) = exp(-norm(xiTheta-cij(5,j))^2/(2*(bij(5,j))^2));
        wHatTheta(j) = WHAT_OLDtheta(j)+...
        Ts*(1/gammaTheta)*bPrime(5)*LandaPrime*abs(eDot(5))^(LandaPrime-1)...
        *sTheta*hTheta(j);
        
        hPsi(j) = exp(-norm(xiPsi-cij(6,j))^2/(2*(bij(6,j))^2));
        wHatPsi(j) = WHAT_OLDpsi(j)+...
        Ts*(1/gammaPsi)*bPrime(6)*LandaPrime*abs(eDot(end))^(LandaPrime-1)*sPsi*hPsi(j);
        
    end
    
    FaultHat(1) = wHatX'*hx;
    FaultHat(2) = wHatY'*hy;
    FaultHat(3) = wHatZ'*hz;
    FaultHat(4) = wHatPhi'*hPhi;
    FaultHat(5) = wHatTheta'*hTheta;
    FaultHat(6) = wHatPsi'*hPsi;
    
end
