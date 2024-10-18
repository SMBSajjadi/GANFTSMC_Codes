function PlotComparedResults_GANFTSMC_NFTSMC()    
    %{
            Comparing the Results of GA-NFTSMC with NFTSMC for both
                                        trajectories
    %}
    
    set(groot, 'DefaultAxesFontName', 'Times New Roman'); % Font name for x and y ticks
    set(groot, 'DefaultAxesFontSize', 12); % Font size for x and y ticks
    set(groot, 'DefaultAxesFontWeight', 'bold'); % Font weight for x and y ticks
    set(groot, 'DefaultLegendFontName', 'Times New Roman'); % Font name for legend
    set(groot, 'DefaultLegendFontSize', 12); % Font size for legend
    set(groot, 'DefaultLegendFontWeight', 'bold'); % Font weight for legend

    t1 = 0:0.001:50;
    Xd = setDesiredTrajectory(t1,1,6);

    t2 = 0:0.001:70;
    Xd2 = setDesiredTrajectory(t2,3,6);

    %% Trajectory 1. Rectangular Trajectory
    
    xGANFTSMC_Rect = load('xGANFTSMC');
    xGANFTSMC_Rect = xGANFTSMC_Rect.x;
    
    uGANFTSMC_Rect = load('uGANFTSMC');
    uGANFTSMC_Rect = uGANFTSMC_Rect.u;
    
    wGANFTSMC = load('wGANFTSMC');
    wGANFTSMC = wGANFTSMC.W;
    
    
    xNFTSMC_Rect = load('xNFTSMC');
    xNFTSMC_Rect = xNFTSMC_Rect.x;
    
    uNFTSMC_Rect = load('uNFTSMC');
    uNFTSMC_Rect = uNFTSMC_Rect.u;
    
    wNFTSMC = load('wNFTSMC');
    wNFTSMC = wNFTSMC.W;
    
    %% Trajectory 2. Combined Trajectory
    
    xGANFTSMC_Com = load('xGANFTSMC_2');
    xGANFTSMC_Com = xGANFTSMC_Com.x;
    
    uGANFTSMC_Com = load('uGANFTSMC_2');
    uGANFTSMC_Com = uGANFTSMC_Com.u;
    
    wGANFTSMC_Com = load('wGANFTSMC_2');
    wGANFTSMC_Com = wGANFTSMC_Com.W;
    
    
    xNFTSMC_Com = load('xNFTSMC_2');
    xNFTSMC_Com = xNFTSMC_Com.x;
    
    uNFTSMC_Com = load('uNFTSMC_2');
    uNFTSMC_Com = uNFTSMC_Com.u;
    
    wNFTSMC_Com = load('wNFTSMC_2');
    wNFTSMC_Com = wNFTSMC_Com.W;

    %% Plot- Rectangular Trajectory

    f612 = figure(612);
    m = 1;
    g = 9.81;
    n = 6;
    
    NameReference = {'x_{d}','y_{d}','z_{d}'};
    YLABELNAME = {'x(m)','y(m)','z(m)'};
                 
    for i=1:n/2
        
        subplot(3,1,i)
        plot(t1,(xGANFTSMC_Rect(i,:)),'b','LineWidth',1.5)
        hold on
        plot(t1,xNFTSMC_Rect(i,:),'m','LineWidth',1)
        hold on
        plot(t1,(Xd(i,:)),'color','r','LineWidth',1,'LineStyle','--')
        grid on
        xlabel('Time (s)','InterPreter','Latex')
        legend('GANFTSMC','NFTSMC',NameReference{i},'FontWeight','bold')
        ylabel(YLABELNAME{i})
%         ylim([0 3])

    end


    f613 = figure(613);
    U_Name = {'F_{T} (N)','\tau_{\phi} (N.m)','\tau_{\theta} (N.m)','\tau_{\psi} (N.m)'};
    uMainZGA_Rect = m*sqrt(uGANFTSMC_Rect(1,:).^2+uGANFTSMC_Rect(2,:).^2+(g+uGANFTSMC_Rect(3,:)).^2);
    uMainZGA_Rect = [uMainZGA_Rect
                                  uGANFTSMC_Rect(4,:)
                                  uGANFTSMC_Rect(5,:)
                                  uGANFTSMC_Rect(6,:)];

    uMainZ_Rect = m*sqrt(uNFTSMC_Rect(1,:).^2+uNFTSMC_Rect(2,:).^2+(g+uNFTSMC_Rect(3,:)).^2);
    uMainZ_Rect = [uMainZ_Rect
                          uNFTSMC_Rect(4,:)
                          uNFTSMC_Rect(5,:)
                          uNFTSMC_Rect(6,:)];

    for i=1:4
    
        subplot(2,2,i)
        plot(t1,uMainZGA_Rect(i,:),'b','LineWidth',1.6)
        grid on
        hold on
        plot(t1,uMainZ_Rect(i,:),'r--','LineWidth',1)
        xlabel('Time(s)','InterPreter','Latex')
        legend('GANFTSMC','NFTSMC','FontWeight','bold')
%         xlim([0 50])
        ylabel(U_Name{i})

        % Set if required
%     if(i==1)
%             ylim([5 40])
%             yticks([5 10 20 30 40])
%         elseif(i==2)
%             ylim([-0.2 0.05])
%         elseif(i==3)
%             ylim([-0.2 0.2])
%             yticks([-0.06 -0.03 0 0.03 0.06])
%         else
%             ylim([0 0.15])
%     end

    end


    f614 = figure(614);
    Name_w = {'\omega_{1} (Rad/s)','\omega_{2} (Rad/s)','\omega_{3} (Rad/s)',...
                        '\omega_{4} (Rad/s)'};
    
    for i=1:4
            
        subplot(2,2,i);
        plot(t1,(wGANFTSMC(i,:)),'b','LineWidth',1.6)
        hold on
        grid on
        plot(t1,(wNFTSMC(i,:)),'r--','LineWidth',1)
        xlabel('Time(s)','InterPreter','Latex')
        legend('GANFTSMC','NFTSMC','FontWeight','bold')
        xlim([0 50])
        ylabel(Name_w{i})

    end























    %% Plot- Combined Trajectory

    f615 = figure(615);
    m = 1;
    g = 9.81;
    n = 6;
    
    NameReference = {'x_{d}','y_{d}','z_{d}'};
    YLABELNAME = {'x(m)','y(m)','z(m)'};
                 
    for i=1:n/2
        
        subplot(3,1,i)
        plot(t2,(xGANFTSMC_Com(i,:)),'b','LineWidth',1.6)
        hold on
        plot(t2,xNFTSMC_Com(i,:),'m','LineWidth',1.1)
        hold on
        plot(t2,(Xd2(i,:)),'r--','LineWidth',1.1,'LineStyle','--')
        grid on
        xlabel('Time (s)','InterPreter','Latex')
        legend('GANFTSMC','NFTSMC',NameReference{i},'FontWeight','bold')
        ylabel(YLABELNAME{i})
%         xlim([0 70])

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


    f616 = figure(616);
    U_Name = {'F_{T} (N)','\tau_{\phi} (N.m)','\tau_{\theta} (N.m)','\tau_{\psi} (N.m)'};
    uMainZGA_Combined = m*sqrt(uGANFTSMC_Com(1,:).^2+uGANFTSMC_Com(2,:).^2+(g+uGANFTSMC_Com(3,:)).^2);
    uMainZGA_Com = [uMainZGA_Combined
                                  uGANFTSMC_Com(4,:)
                                  uGANFTSMC_Com(5,:)
                                  uGANFTSMC_Com(6,:)];

    uMainZ_Combined = m*sqrt(uNFTSMC_Com(1,:).^2+uNFTSMC_Com(2,:).^2+(g+uNFTSMC_Com(3,:)).^2);
    uMainZ_Com = [uMainZ_Combined
                          uNFTSMC_Com(4,:)
                          uNFTSMC_Com(5,:)
                          uNFTSMC_Com(6,:)];

    for i=1:4
    
        subplot(2,2,i)
        plot(t2,uMainZGA_Com(i,:),'b','LineWidth',1.5)
        grid on
        hold on
        plot(t2,uMainZ_Com(i,:),'r--','LineWidth',1.2)
        xlabel('Time(s)','InterPreter','Latex')
        legend('GANFTSMC','NFTSMC','FontWeight','bold')
%         xlim([0 70])
%         xticks([0 10 20 30 40 50 60 70])
        ylabel(U_Name{i})

    end


    f617 = figure(617);
    Name_w = {'\omega_{1} (Rad/s)','\omega_{2} (Rad/s)','\omega_{3} (Rad/s)',...
                        '\omega_{4} (Rad/s)'};
    
    for i=1:4
            
        subplot(2,2,i);
        plot(t2,(wGANFTSMC_Com(i,:)),'b','LineWidth',1.6)
        hold on
        grid on
        plot(t2,(wNFTSMC_Com(i,:)),'r--','LineWidth',1)
        xlabel('Time(s)','InterPreter','Latex')
        legend('GANFTSMC','NFTSMC','FontWeight','bold')
        xlim([0 70])
        ylabel(Name_w{i})
%         ylim([0 600])
%         xticks([0 10 20 30 40 50 60 70])

    end

    f618 = figure(618);
    b = 5.42e-5;    % Drag Force Coefficient
    ylabelF = {'$F_1 (N)$','$F_2 (N)$','$F_3 (N)$','$F_4 (N)$'};
    legendF = {'$F_1$','$F_2$','$F_3$','$F_4$'};

    for i=1:4
        
        subplot(2,2,i)
        plot(t2,b*wGANFTSMC_Com(i,:).^2,'b','LineWidth',1.2)
        grid on
        xlabel('Time (s)','InterPreter','Latex')
        ylabel(ylabelF{i},'InterPreter','Latex')
        legend(legendF{i},'InterPreter','Latex')
%         xlim([0 70])
%         xticks([0 10 20 30 40 50 60 70])

%         if(i==1)
%             ylim([0 5])
%         end
% 
%         if(i==2)
%             ylim([0 3])
%         end
%         if(i==3)
%             ylim([0 5])
%         end
%         if(i==4)
%             ylim([0 3])
%         end

    end

    movegui(f612,'east')
    movegui(f613,'west')
    movegui(f614,'north')
    movegui(f615,'south')
    movegui(f616,'northeast')
    movegui(f617,'northwest')
    movegui(f618,'southwest')

end
