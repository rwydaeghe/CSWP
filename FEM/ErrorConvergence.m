function ErrorConvergence()
    max_modes=5;
    exactTM=[5.7832,15.6528,30.4713,40.3409,45.2616]; exactTE=[24.5515,54.1603,59.0882,88.6971,103.5084];
    TMeigs=zeros(max_modes,1); TEeigs=zeros(max_modes,1); hasConvergedTM=zeros(max_modes,2); hasConvergedTE=zeros(max_modes,2);
    N=0; TM_xaxis=[]; TE_xaxis=[];
    while ~(all(hasConvergedTM(:,1)) & all(hasConvergedTE(:,2)))
        if N<=30
            N=N+1;
        elseif N>30 & N <= 50
            N=N+3; disp(['N = ' num2str(N)])
        elseif N>50 & N <= 70
            N=N+10; disp(['N = ' num2str(N)])
        elseif N>70 & N <= 100
            N=N+15; disp(['N = ' num2str(N)])
        else
            disp('Necessary N for convergence exceeds 100. Aborting procedure.')
            return
        end        
        [TM,TE]=ComputeEigs(N); 
        if size(TM,1)>max_modes
            TM=TM(1:max_modes);
        elseif size(TM,1)<max_modes
            TM=[TM;zeros(max_modes-size(TM,1),1)];
        end
        if size(TE,1)>max_modes
            TE=TE(1:max_modes);
        elseif size(TE,1)<max_modes
            TE=[TE;zeros(max_modes-size(TE,1),1)];
        end
        if ~all(hasConvergedTM(:,1))
            TMeigs(:,N)=TM;
            TM_xaxis=[TM_xaxis,N]; 
        end
        if ~all(hasConvergedTE(:,1))
            TEeigs(:,N)=TE;
            TE_xaxis=[TE_xaxis,N];
        end
        for mode=1:size(TM,1)
            if abs(TM(mode)-exactTM(mode))/TM(mode)<0.01 & hasConvergedTM(mode)==0
                hasConvergedTM(mode,:)=[N,TM(mode)];
            end
        end
        for mode=1:size(TE,1)
            if abs(TE(mode)-exactTE(mode))/TE(mode)<0.01 & hasConvergedTE(mode)==0
                hasConvergedTE(mode,:)=[N,TE(mode)];
            end
        end
    end
    legendLabels={};   
    for mode=1:max_modes
        legendLabels{end+1}=['mode #' num2str(mode)];
    end
    figure('Name','TM resonant frequencies k^2')
    hold on  
    plot(TM_xaxis,TMeigs(:,TM_xaxis),'-x')
    scatter(hasConvergedTM(:,1),hasConvergedTM(:,2))    
    title('Error convergence of resonant frequencies as function of mesh size $N$','interpreter','latex')
    legend(legendLabels,'FontSize',8,'Location','best')
    xlabel('$N$','interpreter','latex')
    ylabel('$k^2$','interpreter','latex')
    hold off
    
    figure('Name','TE resonant frequencies k^2')    
    hold on
    plot(TE_xaxis,TEeigs(:,TE_xaxis),'-x')
    scatter(hasConvergedTE(:,1),hasConvergedTE(:,2))        
    title('Error convergence of resonant frequencies as function of mesh size $N$','interpreter','latex')
    legend(legendLabels,'FontSize',8,'Location','best')
    xlabel('$N$','interpreter','latex')
    ylabel('$k^2$','interpreter','latex')    
    hold off
    
    hasConvergedTM
    hasConvergedTE
end
