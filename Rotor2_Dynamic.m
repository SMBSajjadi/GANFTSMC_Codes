function [xDot, W, Disturbance,wStar, FvectorFault] = Rotor2_Dynamic(t,x,u,FAULT_ANGLES)
    %% A. System Parameters
    
        g = 9.81;         % Gravity Acceleration
        L = 0.47/2;      % One-Half Length. Full Length equals 2*L = 47cm
        m = 1;            % Mass of the Quadrotor
        
        Ix = 0.0081;        % X Axis Moment of Intertia
        Iy = Ix;               % Y Axis Moment of Intertia
        Iz = 0.0142;        % Z Axis MOMENT of Intertia
        JTP = 10.4e-5;    
        
        b = 5.42e-5;    % Drag Force Coefficient
        d = 1.1e-6;      % Drag Torque Coefficient
        
        Kf = 1e-6;
        Kt = 1.2e-6;
        
        AlphaAngle = FAULT_ANGLES(1);
        BetaAngle = FAULT_ANGLES(2);
        GammaAngle = FAULT_ANGLES(3);
        
        %% Control Part
    
        Transform_Matrix = [b       b       b       b
                                        0    -b*L    0       b*L
                                       -b*L   0      b*L    0
                                        d     -d      d       -d];

        uThrust = m*sqrt(u(1)^2+u(2)^2+(g+u(3))^2);
        Sol_Vector = [uThrust u(4) u(5) u(6)]';        % [Uz, Uphi, Utheta, Upsi]'
        w2_2 = uThrust/(4*b) - u(6)/(4*d) - u(4)/(2*b*L);
        Squared_W = (linsolve(Transform_Matrix,Sol_Vector));
    
        w1s = Squared_W(1);
        w2s = Squared_W(2);
        w3s = Squared_W(3);
        w4s = Squared_W(4);
        
        W = real(sqrt([w1s w2s w3s w4s]));
        w1 = W(1);
        w2 = W(2);
        w3 = W(3);
        w4 = W(4);
        
    %% B. State Vector
    
%     X = x(1);
%     y = x(2);
%     z = x(3);
    
    phi = x(4);
    theta = x(5);
    say = x(6);
    
    xDot = x(7);
    yDot = x(8);
    zDot = x(9);
    
    phiDot = x(10);
    thetaDot = x(11);
    sayDot = x(12);
    
    c = @(x) cos(x);                    % Cosine Function
    s = @(x) sin(x);                     % Sinusoidal Function
    
    wStar = (w1 + w3 - w2 - w4);        % Disturbance
     
    %% C. Fault Injection 

        f1 = s(AlphaAngle)*s(GammaAngle);
        f2 = -c(GammaAngle)*s(BetaAngle) + s(GammaAngle)*c(BetaAngle)*c(AlphaAngle);
        f3 = c(BetaAngle)*c(GammaAngle) + c(AlphaAngle)*s(BetaAngle)*s(GammaAngle)-1;
        f4 = f2*s(BetaAngle)-(1+f3)*c(BetaAngle)+1;
        f5 = f1*s(BetaAngle);
        f6 = f1*c(BetaAngle);
        
        FvectorFault = [f1 f2 f3 f4 f4 f6]';
            
        ufx = (b/m)*w2_2*(+f1*(c(theta)*c(say)) + f2*(c(say)*s(phi)*s(theta) - c(phi)*s(say)) + f3*(s(phi)*s(say) + c(phi)*c(say)*s(theta)));
        ufy = (b/m)*w2_2*(+f1*(c(theta)*s(say)) + f2*(c(phi)*c(say) + s(phi)*s(theta)*s(say)) + f3*(-c(say)*s(phi) + c(phi)*s(say)*s(theta)));
        ufz = (b/m)*w2_2*(-f1*s(theta) + f2*c(theta)*s(phi) + f3*c(phi)*c(theta));

        ufPhi = (JTP*w2/Ix)*(sayDot*f2 - f3*thetaDot) + (1/Ix)*w2_2*(b*L*f4 + f1*d);
        ufTheta = (JTP/Iy)*w2*(sayDot*f1 + f3*phiDot) + (1/Iy)*(w2_2)*(-f2*d + L*b*f5);
        ufSay = (JTP*w2)/Iz*(-f1*thetaDot - f2*phiDot) + (1/Iz)*(w2_2)*(-f3*d - f6*L*b);
        
        %% Unkonown Disturbance to be Estimated Using RNN
        
        Disturbance = [ufx      ufy         ufz         ufPhi       ufTheta         ufSay]';    
       
        %% State Space 
                     
           xDoubleDot = u(1)-Kf*xDot/m+ufx;
           yDoubleDot = u(2)-Kf*yDot/m+ufy;
           zDoubleDot = u(3)-Kf*zDot/m+ufz;
           
           phiDoubleDot = ((Iy-Iz)/Ix)*thetaDot*sayDot+JTP*thetaDot*wStar/Ix+u(4)/Ix-Kt*L*phiDot/Ix+ufPhi;
           thetaDoubleDot = ((Iz-Ix)/Iy)*phiDot*sayDot-JTP*phiDot*wStar/Iy+u(5)/Iy-(Kt*L/Iy)*thetaDot+ufTheta;
           psiDoubleDot = ((Ix-Iy)/Iz)*phiDot*thetaDot+u(6)/Iz-(Kt*L/Iz)*sayDot+ufSay;
           
              xDot =  [x(7:12)
                           xDoubleDot
                           yDoubleDot
                           zDoubleDot
                           phiDoubleDot
                           thetaDoubleDot
                           psiDoubleDot];
 
end
