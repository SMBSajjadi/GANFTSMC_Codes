function [u,S,e,eDot] = ...
    RBFNNTSMC_GENETIC(wStar,x,n,XD,XDotD,XDoubleDotD, b,bPrime, Landa,LandaPrime,...
                                      K,FHAT_ST)

    FHAT_X = FHAT_ST(1);
    FHAT_Y = FHAT_ST(2);
    FHAT_Z = FHAT_ST(3);
    FHAT_PHI = FHAT_ST(4);
    FHAT_THETA = FHAT_ST(5);
    FHAT_PSI = FHAT_ST(6);

    xDot = x(7);
    yDot = x(8);
    zDot =x(9);
    phiDot = x(10);
    thetaDot = x(11);
    psiDot = x(12);
    
    L = 0.47/2;          % One-Half Length. Full Length equals 2*L = 47cm
    m = 1;                % Mass of the Quadrotor

    Ix = 0.0081;        % X Axis Moment of Intertia
    Iy = Ix;               % Y Axis Moment of Intertia
    Iz = 0.0142;        % Z Axis MOMENT of Intertia
    JTP = 10.4e-5;    

    Kf = 1e-6;
    Kt = 1.2e-6;
    
    e = x(1:n) - XD(1:n);
    eDot = x(n+1:end) - XDotD(1:n);
    ex = e(1);
    ey = e(2);
    ez = e(3);
    ePhi = e(4);
    eTheta = e(5);
    ePsi = e(6);
    
    eDotx = eDot(1);
    eDoty = eDot(2);
    eDotz = eDot(3);
    eDotPhi = eDot(4);
    eDotTheta = eDot(5);
    eDotPsi = eDot(6);
    
    %% Nonsingular Sliding Surfaces
    
    Sx = ex +...
        b(1)*((abs(ex))^Landa)*sign(ex) +...
        bPrime(1)*(abs(eDotx)^LandaPrime)*sign(eDotx);
    
    Sy = ey +...
        b(2)*((abs(ey))^Landa)*sign(ey) +...
        bPrime(2)*(abs(eDoty)^LandaPrime)*sign(eDoty);
    
    Sz = ez +...
        b(3)*((abs(ez))^Landa)*sign(ez) +...
        bPrime(3)*(abs(eDotz)^LandaPrime)*sign(eDotz);
    
    Sphi = ePhi +...
        b(4)*((abs(ePhi))^Landa)*sign(ePhi) +...
        bPrime(4)*(abs(eDotPhi)^LandaPrime)*sign(eDotPhi);
    
    Stheta = eTheta +...
        b(5)*((abs(eTheta))^Landa)*sign(eTheta) +...
        bPrime(5)*(abs(eDotTheta)^LandaPrime)*sign(eDotTheta);
    
    Spsi = ePsi +...
              b(6)*((abs(ePsi))^Landa)*sign(ePsi) +...
              bPrime(6)*(abs(eDotPsi)^LandaPrime)*sign(eDotPsi);
   
    %% X-Direction
    
    RegulationPartX = -bPrime(1)^-1*LandaPrime^-1*(abs(eDotx)^(2-LandaPrime))*sign(eDotx)*...
        (1+b(1)*Landa*abs(ex)^(Landa-1));
    TrackingPartX = Kf*xDot/m+XDoubleDotD(1)-K(1)*Sx-K(2)*tanh(Sx)-FHAT_X;
    ux = (RegulationPartX+TrackingPartX);
    
    %% Y-Direction
    
    RegulationPartY = -bPrime(2)^-1*LandaPrime^-1*(abs(eDoty)^(2-LandaPrime))*sign(eDoty)*...
    (1+b(2)*Landa*abs(ey)^(Landa-1));
    TrackingPartY = Kf*yDot/m+XDoubleDotD(2)-K(1)*Sy-K(2)*tanh(Sy)-FHAT_Y;
    uy = (RegulationPartY+TrackingPartY);
    
    %% Z-Direction
    
    RegulationPartZ = -bPrime(3)^-1*LandaPrime^-1*(abs(eDotz)^(2-LandaPrime))*sign(eDotz)*...
    (1+b(3)*Landa*abs(ez)^(Landa-1));
    TrackingPartZ = Kf*zDot/m+XDoubleDotD(3)-K(1)*Sz-K(2)*tanh(Sz)-FHAT_Z;
    uz = (RegulationPartZ+TrackingPartZ);
    
    %% PHI-Direction
    
    RegulationPartPhi = bPrime(4)^-1*LandaPrime^-1*(abs(eDotPhi)^(2-LandaPrime))*sign(eDotPhi)*...
    (1+b(4)*Landa*abs(ePhi)^(Landa-1));
    TrackingPartPhi = -XDoubleDotD(4)+((Iy-Iz)/Ix)*thetaDot*psiDot+...
        JTP*thetaDot*wStar/Ix-Kt*L*phiDot/Ix+K(1)*Sphi+K(2)*tanh(Sphi)+FHAT_PHI;
    uPhi = -Ix*(RegulationPartPhi+TrackingPartPhi);
    
    %% Theta-Direction
    
    RegulationPartTheta = bPrime(5)^-1*LandaPrime^-1*(abs(eDotTheta)^(2-LandaPrime))*sign(eDotTheta)*...
    (1+b(5)*Landa*abs(eTheta)^(Landa-1));
    TrackingPartTheta = -XDoubleDotD(5)+((Iz-Ix)/Iy)*phiDot*psiDot...
        -JTP*phiDot*wStar/Iy-Kt*L*thetaDot/Iy+K(1)*Stheta+K(2)*tanh(Stheta)+FHAT_THETA;
    uTheta = -Iy*(RegulationPartTheta+TrackingPartTheta);
    
    %% PSI-Direction
    
    RegulationPartPsi = bPrime(6)^-1*LandaPrime^-1*(abs(eDotPsi)^(2-LandaPrime))*sign(eDotPsi)*...
    (1+b(6)*Landa*abs(ePsi)^(Landa-1));
    TrackingPartPsi = -XDoubleDotD(6)+((Ix-Iy)/Iz)*phiDot*thetaDot-Kt*L*psiDot/Iz...
                              +K(1)*Spsi+K(2)*tanh(Spsi)+FHAT_PSI;
    uPsi = -Iz*(RegulationPartPsi+TrackingPartPsi);

    %% NFTSMC Cont.
    
    S = ([Sx, Sy, Sz, Sphi, Stheta, Spsi]');
    u = ([ux, uy, uz, uPhi, uTheta, uPsi]');
    
end
