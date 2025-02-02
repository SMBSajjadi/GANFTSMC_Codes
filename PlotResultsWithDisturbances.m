function PlotResultsWithDisturbances()

    clc;clear;close all;
    set(groot, 'DefaultAxesFontName', 'Times New Roman'); % Font name for x and y ticks
    set(groot, 'DefaultAxesFontSize', 12); % Font size for x and y ticks
    set(groot, 'DefaultAxesFontWeight', 'bold'); % Font weight for x and y ticks
    set(groot, 'DefaultLegendFontName', 'Times New Roman'); % Font name for legend
    set(groot, 'DefaultLegendFontSize', 12); % Font size for legend
    set(groot, 'DefaultLegendFontWeight', 'bold'); % Font weight for legend

    t = 0:0.001:50;
    xGANFTSMC = load('xGANFTSMC.mat');
    xGANFTSMC = xGANFTSMC.x;

    xDOBTSMC = load('xDOBTSMC.mat');
    xDOBTSMC = xDOBTSMC.x;

    XD = load('XD','XD');
    XD = XD.XD;

    f1313 = figure(1);
    NamesPosition = {'Method [33]','Proposed'};
    XDNamesPosition = {'x_{d}','y_{d}','z_{d}'};
    yLabelPos = {'x (m)','y (m)','z (m)'};

    for i=1:3
        
        subplot(3,1,i)
        plot(t,xDOBTSMC(i,:),'color',[0.5 0.1 0.7],'LineWidth',1.3)
        hold on
        grid on
        plot(t,xGANFTSMC(i,:),'r','LineWidth',1.1)
        hold on
        plot(t,XD(i,:),'k--','LineWidth',1)
        legend(NamesPosition{1},NamesPosition{2},XDNamesPosition{i})
        xlabel('Time (s)','InterPreter','Latex')
        ylabel(yLabelPos{i})

        if(i==1)
            ylim([-2 6])
        elseif(i==2)
            ylim([-2 6])
        else
            ylim([-4 4])
        end

    end


    f1314 = figure(2);
    XDNamesOri = {'\phi_{d}','\theta_{d}','\psi_{d}'};
    j = 0;
    yLabelNamesOri = {'\phi (Degree)','\theta (Degree)','\psi (Degree)'};

    for i=4:6
        
        j = j + 1;
        subplot(3,1,j)
        plot(t,rad2deg(xDOBTSMC(i,:)),'color',[0.5 0.1 0.7],'LineWidth',1.3)
        hold on
        grid on
        plot(t,rad2deg(xGANFTSMC(i,:)),'r','LineWidth',1.1)
        hold on
        plot(t,rad2deg(XD(i,:)),'k--','LineWidth',1)
        legend(NamesPosition{1},NamesPosition{2},XDNamesOri{j})
        ylabel(yLabelNamesOri{j})
        xlabel('Time (s)','InterPreter','Latex')

        if(i==6)
            ylim((180/pi)*[-0.05 0.1])
        elseif(i==5)
            ylim((180/pi)*[-0.6 0.2])
        else
            ylim((180/pi)*[-0.2 0.4])
        end

        
    end

    f1315 = figure(3);
    plot3(xDOBTSMC(1,:),xDOBTSMC(2,:),xDOBTSMC(3,:),'color',[0.5 0.1 0.7],'LineWidth',1.3)
    hold on
    grid on
    plot3(xGANFTSMC(1,:),xGANFTSMC(2,:),xGANFTSMC(3,:),'r','LineWidth',1.1)
    hold on
    grid on
    plot3(XD(1,:),XD(2,:),XD(3,:),'k--','LineWidth',1)
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    legend('Method [33]','Proposed','Reference')

    movegui(f1313,'east')
    movegui(f1314,'west')
    movegui(f1315,'north')

end