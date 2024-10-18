function plotResults(t,x,Xd,XDotD,u,n,Omega,S)
    %% A. Plot the States

    set(groot, 'DefaultAxesFontName', 'Times New Roman'); % Font name for x and y ticks
    set(groot, 'DefaultAxesFontSize', 12); % Font size for x and y ticks
    set(groot, 'DefaultAxesFontWeight', 'bold'); % Font weight for x and y ticks
    set(groot, 'DefaultLegendFontName', 'Times New Roman'); % Font name for legend
    set(groot, 'DefaultLegendFontSize', 12); % Font size for legend
    set(groot, 'DefaultLegendFontWeight', 'bold'); % Font weight for legend

    f1 = figure(1);
    m = 1;
    g = 9.81;
    
    Name = {'x','y','z'};
    NameReference = {'x_{d}','y_{d}','z_{d}'};
    YLABELNAME = {'x(m)','y(m)','z(m)'};
                 
    for i=1:n/2
        
        subplot(3,1,i)
        plot(t,(x(i,:)),'b','LineWidth',1.5)
        hold on
        plot(t,(Xd(i,:)),'color','r','LineWidth',1.7,'LineStyle','--')
        grid on
        xlabel('Time (s)','InterPreter','Latex')
        legend(Name{i},NameReference{i},'FontWeight','bold')
        ylabel(YLABELNAME{i})
%         ylim([0 3])
%         if(i==1)
%             ylim([-1 4])
%         end
%         if(i==2)
%             ylim([-1 3])
%         end
%         if(i==3)
%             ylim([0 3])
%         end

    end
    
    f8 = figure(8);
    Name = {'\phi','\theta','\psi'};
    NameReference = {'\phi_{d}','\theta_{d}','\psi_{d}'};
    ii = 0;
    
    for i=n/2+1:n
        
        ii = ii + 1;
        subplot(3,1,ii)
        plot(t,180/pi*((x(i,:))),'b','LineWidth',1.3)
        hold on
        plot(t,180/pi*((Xd(i,:))),'color','r','LineWidth',1.1,'LineStyle','-.')
        grid on
        xlabel('Time (s)','InterPreter','Latex')
        legend(Name{ii},NameReference{ii})
        ylabel([Name{ii}, num2str(' (Degree)')])
        
        if(ii==1)
            ylim([-40 40])
        end

        if(ii==2)
            ylim([-40 40])
        end

        if(ii==3)
            ylim([0 40])
            yticks([0 10 20 30 40])
        end
        
    end

    %% B. Plot The Remaining States
    
%     f2 = figure(2);
%     Name2 = {'$\dot{x}$','$\dot{y}$','$\dot{z}$',...
%              '$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$'};
%          
%     Name2Reference = {'$\dot{x}_{d}$','$\dot{y}_{d}$','$\dot{z}_{d}$',...
%                                    '$\dot{\phi}_{d}$','$\dot{\theta}_{d}$','$\dot{\psi}_{d}$'};
%     yLabelName2 = {'$\dot{x} (m/s)$','$\dot{y} (m/s)$','$\dot{z} (m/s)$',...
%                              '$\dot{\phi} (rad/s)$','$\dot{\theta} (rad/s)$','$\dot{\psi} (rad/s)$'};
%                  
%     iter = 0;
%          
%     for i=n+1:2*n
%         
%         iter = iter + 1;
%         subplot(3,2,iter)
%         plot(t,x(i,:),'LineWidth',2)
%         grid on
%         hold on
%         plot(t,XDotD(i-n,:),'color','k','LineWidth',1.75,'LineStyle','--')
%         xlabel('Time (s)','InterPreter','Latex')
%         legend(Name2{iter},Name2Reference{iter},'InterPreter','latex')
%         ylabel(yLabelName2{iter},'InterPreter','Latex')
% 
%     end
    
    %% Control Signals Plot
    
    f3 = figure(3);
    U_Name = {'F_{T}','\tau_{\phi}','\tau_{\theta}','\tau_{\psi}'};
    uMainZ = m*sqrt(u(1,:).^2+u(2,:).^2+(g+u(3,:)).^2);
    uMain = [uMainZ
                  u(4,:)
                  u(5,:)
                  u(6,:)];
    Ulabel = {'Total Thrust (N)','Roll Torque (N.m)','Pitch Torque (N.m)','Yaw Torque (N.m)'};

    for i=1:4
    
        subplot(2,2,i)
        plot(t,uMain(i,:),'r','LineWidth',1.5)
        xlabel('Time(s)','InterPreter','Latex')
        grid on
        legend(U_Name{i},'FontWeight','bold')
        xlim([0 50])
        ylabel(Ulabel{i})
        
        % Trajectory 2

%         if(i==1)
%             ylim([6 16])
%         end
%         if(i==2)
%             ylim([-0.2 0.05])
%         end
% 
%         if(i==3)
%             ylim([-0.1 0.1])
%         end
%         if(i==4)
%             yticks([0 0.02 0.05 0.08 0.1])
%         end


        % Trajectory 1
        
%         if(i==1)
%             ylim([5 20])
%         elseif(i==2)
%             ylim([-0.2 0.05])
%         elseif(i==3)
%             ylim([-0.06 0.06])
%             yticks([-0.06 -0.03 0 0.03 0.06])
%         else
%             ylim([0 0.15])
%         end
        
    end
    
%     f9 = figure(9);
%     ColorU = {'r','b',[0.5 0.1 0.3]};
%     U_Name = {'\tau_{\phi}','\tau_{\theta}','\tau_{\psi}'};
%     
%     subIter = 0;
%     
%     for i=size(u,1)/2+1:n
%     
%         subIter = subIter + 1;
%         subplot(2,2,subIter)
%         if(i==n)
%             subplot(2,2,[subIter,subIter+1])
%         end
%         
%         plot(t,u(i,:),'LineWidth',2,'Color',ColorU{subIter})
%         xlabel('Time(s)','InterPreter','Latex')
%         grid on
%         legend(U_Name{subIter})
%         
%     end
   
    f4 = figure(4);
    Name_w = {'\omega_{1}','\omega_{2}','\omega_{3}','\omega_{4}'};
    
    for i=1:size(Omega,1)
            
        subplot(2,2,i);
        plot(t,real(Omega(i,:)),'LineWidth',1.6)
        xlabel('Time(s)','InterPreter','Latex')
        grid on
        legend(Name_w{i},'FontWeight','bold')
        xlim([0 50])
        ylabel('Angular Velocity (rad/s)')
        
%         if(i==1 || i==3)
%             ylim([100 350])
%         end
% 
%         if(i==4)
%                 ylim([100 300])
%         end
    end
   
    %% Desired 3D Trajectory
    
    f5 = figure(5);
    X = x(1,:);
    y = x(2,:);
    z = x(3,:);
    
    plot3(X,y,z,'LineWidth',1.1,'Color','r')
    xlabel('x (m)','InterPreter','Latex')
    ylabel('y (m)','InterPreter','Latex')
    zlabel('z (m)','InterPreter','Latex')
    hold on
    plot3(Xd(1,:),Xd(2,:),Xd(3,:),'Color','b','LineWidth',1.3,'LineStyle','--')
    grid on
    legend('Flight Trajectory','Target Trajectory','InterPreter','Latex')
    
    %% Tracking Error
    
%     f6 = figure(6);
%     Name_e = {'e_{x}','e_{y}','e_{z}',...
%                       'e_{\phi}','e_{\theta}','e_{\psi}'};
%     for i=1:n
%         
%         subplot(3,2,i)
%         RandomColor = unifrnd(0,1,1,3);
%         plot(t,x(i,:)-Xd(i,:),'LineWidth',2,'Color',RandomColor);
%         xlabel('Time (s)','FontWeight','Bold')
%         legend(Name_e{i})
%         grid on
%         
%     end
    
    %% Plot Sliding Surfaces
    
%     f7 = figure(7);
%     j = 0;
%     
%     for i=1:size(S,1)
%        
%         j = j + 1;
%         subplot(3,2,j)
%         RandomColorVector = unifrnd(0,1,1,3);
%         plot(t,(S(i,:)),'color',RandomColorVector,'LineWidth',2)
%         grid on
%         xlabel('Time (s)','InterPreter','Latex')
%         legend(['S_{' num2str(i) '}'])
%         
%     end
    
    %% Move Figures
    
    movegui(f1,'center')
%     movegui(f2,'east')
    movegui(f3,'west')
    movegui(f4,'north')
    movegui(f5,'south')
%     movegui(f6,'northeast')
%     movegui(f7,'southwest')
    movegui(f8,'northwest')
%     movegui(f9,'southwest')
    
end