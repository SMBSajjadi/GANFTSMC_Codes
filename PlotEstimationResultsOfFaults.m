function PlotEstimationResultsOfFaults(t,MainFaults,EstimatedFaults)

    f1 = figure(11);
    n = 6;
    
    FaultName = {'$f_{stx}$','$f_{sty}$','$f_{stz}$','$f_{st\phi}$','$f_{st\theta}$','$f_{st\psi}$'};
    EstimatedFaultName = {'$\hat{f}_{stx}$', '$\hat{f}_{sty}$', '$\hat{f}_{stz}$', ...
                      '$\hat{f}_{st{\phi}}$', '$\hat{f}_{st{\theta}}$', '$\hat{f}_{st{\psi}}$'};

    ErrorName = {'$e_{fst_x}$','$e_{fst_y}$','$e_{fst_z}$',...
                           '$e_{fst_\phi}$','$e_{fst_\theta}$','$e_{fst_\psi}$'};
                  
    EstimationError = MainFaults - EstimatedFaults;
    for i=1:n
        
        subplot(3,2,i)
        plot(t,MainFaults(i,:),'LineWidth',2)
        hold on
        grid on
        plot(t,EstimatedFaults(i,:),'--','LineWidth',2)
        xlabel('Time (s)','InterPreter','Latex')
        legend(FaultName{i},EstimatedFaultName{i},'InterPreter','Latex')

    end
    
    f2 = figure(15);
    yLabelNamesError = {'$e_{fst_x} (m/s^2)$','$e_{fst_y} (m/s^2)$','$e_{fst_z} (m/s^2)$',...
                     '$e_{fst_\phi} (rad/s^2)$','$e_{fst_\theta} (rad/s^2)$','$e_{fst_\psi} (rad/s^2)$'};
    for i=1:n
        
        subplot(3,2,i)
        plot(t,EstimationError(i,:),'color',[0.8 0.2 0.4],'LineWidth',1.2)
        grid on
        xlabel('Time (s)','InterPreter','Latex')
        legend(ErrorName{i},'InterPreter','Latex')
        xlim([0 70])
        ylabel(yLabelNamesError{i},'InterPreter','Latex')

        if(i==1)
            ylim([-0.2 1])
        end

        if(i==2)
            ylim([-0.5 1.5])
        end
        
        if(i==3)
            ylim([-0.5 0.5])
            yticks([-0.5 -0.25 0 0.25 0.5])
        end

        if(i==4)
            yticks([-15 -10 -5 0 5 10 15])
        end
        
        if(i==6)
            ylim([-10 5])
            yticks([-10 -5 0 5])
        end

    
    end

     movegui(f1,'south')    
     movegui(f2,'east')   

end