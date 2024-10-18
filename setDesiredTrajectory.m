function [XD, XDotD, XDoubleDotD] = setDesiredTrajectory(t,CASE,n)

        if(isempty(CASE) || CASE <=0 || CASE > 3)  
            CASE = 1;
        end
        
        if(nargin < 3)
            n = 6;
        end
        
        switch CASE
            case 1
            %% First Trajectory
            
            xd = @(t) (t<=10)*2 + (t>10 & t<30)*1 + (t>=30)*2;
            yd = @(t) (t<=20)*2 + (t>20 & t<=40)*1 + (t>40)*2;
            zd = @(t) (t<=50)*2 + (t>50)*0;
            psid = @(t) (t<=50)*0.5 + (t>50)*0;

            xd = xd(t);
            yd = yd(t);
            zd = zd(t);
            psid = psid(t);

            phid = zeros(size(t));
            thetad = phid;
            
            XD = [xd
                      yd
                      zd
                      phid
                      thetad
                      psid];
                  
            XDotD = zeros(n,numel(t));
            XDoubleDotD = zeros(n,numel(t));
            
            case 2
            %% Second TRAJECTORY  
            
%             xdesired = @(t) (t<=10)*0 + (t>10 & t<30)*0.3*cos(pi*t/6) + (t>=30)*0;
%             ydesired = @(t) (t<=20)*0 + (t>20 & t<35)*0.3*sin(pi*t/6) + (t>=35)*0;
%             zdesired = @(t) (t<=15)*0.8 + (t>15 & t<35)*0.6 + (t>=35)*1;
            
            xdesired = @(t) 0.5*cos(pi*t/20);
            ydesired = @(t) 0.5*sin(pi*t/20);
            zdesired = @(t) 2-0.5*cos(pi*t/20);
            
            xd = zeros(size(t));
            yd = xd;
            zd = xd;
            
            for i=1:numel(t)
                
                xd(i) = xdesired(t(i));
                yd(i) = ydesired(t(i));
                zd(i) = zdesired(t(i));
                
            end

            psid = zeros(size(t));
            phid = zeros(size(t));
            thetad = phid;
 
            xDotd = @(t) -0.5*pi/20*sin(pi*t/20);
            yDotd = @(t) 0.5*pi/20*cos(pi*t/20);
            zDotd = @(t)  +0.5*pi/20*sin(pi*t/20);  
            
            xDotDesired = zeros(size(t));
            yDotDesired = xDotDesired;
            zDotDesired = zeros(size(t));
            
            for k = 1:numel(t)
                
               xDotDesired(k) = xDotd(t(k));
               yDotDesired(k) = yDotd(t(k));
               zDotDesired(k) = zDotd(t(k));
               
            end
            
            phiDotDesired = zeros(size(t));
            thetaDotDesired = zeros(size(t));
            sayDotDesired = zeros(size(t));
            
            xDoubleDotD = @(t) -(0.5*pi/20)^2*cos(pi*t/20);
            yDoubleDotD = @(t) -(0.5*pi/20)^2*sin(pi*t/20);
            zDoubleDotD = @(t) +(0.5*pi/20)^2*cos(pi*t/20);
            
            xDoubleDotDesired = zeros(size(t));
            yDoubleDotDesired = zeros(size(t));
            zDoubleDotDesired = zeros(size(t));
            
            for i=1:numel(t)
                
                xDoubleDotDesired(i) = xDoubleDotD(t(i));
                yDoubleDotDesired(i) = yDoubleDotD(t(i));
                zDoubleDotDesired(i) = zDoubleDotD(t(i));
                
            end
            
            phiDoubleDotDesired = zeros(size(t));
            thetaDoubleDotDesired = zeros(size(t));
            sayDoubleDotDesired = zeros(size(t));
                  
            XD = [xd
                      yd
                      zd
                      phid
                      thetad
                      psid];
            
            XDotD = [xDotDesired
                           yDotDesired
                           zDotDesired
                           phiDotDesired
                           thetaDotDesired
                           sayDotDesired];
                       
            XDoubleDotD = [xDoubleDotDesired
                                     yDoubleDotDesired
                                     zDoubleDotDesired
                                     phiDoubleDotDesired
                                     thetaDoubleDotDesired
                                     sayDoubleDotDesired];
            
            case 3
             %% 3rD Trajectory
                
                 xd = @(t) (t<=4*pi)*(0.5*cos(t/2)) +...
                                (t>=4*pi & t<20)*0.5 + (t>=20 & t<30)*(0.25*t-4.5)+...
                                (t>=30)*3;
                 
                 yd = @(t) (t<=4*pi)*(0.5*sin(t/2)) +...
                     (t>=4*pi & t<20)*(0.25*t-3.14) +...
                     (t>=20 & t<30)*(5-pi)+...
                     (t>=30 & t<40)*(-0.2358*t+8.94) + ...
                     (t>=40)*(-0.5);
                 
                 zd = @(t) (t<4*pi)*(0.125*t+1) + ...
                                (t>=4*pi & t<40)*(0.5*pi+1) +...
                                (t>=40)*exp(-0.2*t+8.944);
                 
                xdesired = zeros(size(t));
                ydesired = xdesired;
                zdesired = xdesired;

            for i=1:numel(t)
                
                xdesired(i) = xd(t(i));
                ydesired(i) = yd(t(i));
                zdesired(i) = zd(t(i));
                
            end
            
            xd = xdesired;
            yd = ydesired;
            zd = zdesired;

            psid = zeros(size(t)); sayDotDesired = psid;    sayDoubleDotDesired = psid;
            phid = zeros(size(t)); phiDotDesired = psid;    phiDoubleDotDesired = psid;
            thetad = phid;         thetaDotDesired = psid;    thetaDoubleDotDesired = psid;
                
            xDotDesired = @(t) (t<=4*pi)*(-0.25*sin(t/2)) +...
                                         (t>=4*pi & t<20)*0 + (t>=20 & t<30)*(0.25)+...
                                         (t>=30)*0;
                   
            yDotDesired = @(t) (t<=4*pi)*(0.25*cos(t/2)) +...
                                         (t>=4*pi & t<20)*(0.25) +...
                                         (t>=20 & t<30)*(0)+...
                                         (t>=30 & t<40)*(-0.2358) + ...
                                         (t>=40)*(0);
                 
            zDotDesired = @(t) (t<4*pi)*(0.125) + ...
                                          (t>=4*pi & t<40)*(0) +...
                                          (t>=40)*(-0.2*exp(-0.2*t+8.944));
                 
                 xDOTDESIRED = zeros(size(t));
                 yDOTDESIRED = zeros(size(t));
                 zDOTDESIRED = zeros(size(t));
                 
             for i=1:numel(t)
                
                xDOTDESIRED(i) = xDotDesired(t(i));
                yDOTDESIRED(i) = yDotDesired(t(i));
                zDOTDESIRED(i) = zDotDesired(t(i));
                
             end
             
              xDoubleDotDesired = @(t) (t<=4*pi)*(-0.125*cos(t/2)) +...
                                                     (t>=4*pi & t<20)*0 + (t>=20 & t<30)*(0)+...
                                                     (t>=30)*0;
                 
                 
             yDoubleDotDesired = @(t) (t<=4*pi)*(-0.125*sin(t/2)) +...
                                                     (t>=4*pi & t<20)*(0) +...
                                                     (t>=20 & t<30)*(0)+...
                                                     (t>=30 & t<40)*(0) + ...
                                                     (t>=40)*(0);
                 
             zDoubleDotDesired = @(t) (t<4*pi)*(0) + ...
                                                    (t>=4*pi & t<40)*(0) +...
                                                    (t>=40)*(0.04*exp(-0.2*t+8.944));
                 
                 xDOUBLEDOTDESIRED = zeros(size(t)); yDOUBLEDOTDESIRED = zeros(size(t));
                 zDOUBLEDOTDESIRED = zeros(size(t));
                 
              for i=1:numel(t)
                
                xDOUBLEDOTDESIRED(i) = xDoubleDotDesired(t(i));
                yDOUBLEDOTDESIRED(i) = yDoubleDotDesired(t(i));
                zDOUBLEDOTDESIRED(i) = zDoubleDotDesired(t(i));
                
              end
            
              XD = [xd
                        yd
                        zd
                        phid
                        thetad
                        psid];
              
              XDotD = [xDOTDESIRED
                             yDOTDESIRED
                             zDOTDESIRED
                             phiDotDesired
                             thetaDotDesired
                             sayDotDesired];
              
               XDoubleDotD = [xDOUBLEDOTDESIRED
                                        yDOUBLEDOTDESIRED
                                        zDOUBLEDOTDESIRED
                                        phiDoubleDotDesired
                                        thetaDoubleDotDesired
                                        sayDoubleDotDesired];
               
        end

end
