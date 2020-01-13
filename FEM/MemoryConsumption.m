function MemoryConsumption(scenario)
    max_modes=2;
    if scenario=="empty"
        exactTM=[5.7832,15.6528,30.4713,40.3409,45.2616]; exactTE=[24.5515,54.1603,59.0882,88.6971,103.5084];
    elseif scenario=="half-filled"
        exactTM=[1.8077,2.9359,3.5500,5.0948,5.2435].^2; exactTE=[3.2996,4.8998,5.6206,6.5751,6.9783].^2;
    elseif scenario=="partly-filled"
        exactTM=[-300,-300,-300,-300,-300].^2; exactTE=[-300,-300,-300,-300,-300].^2; %will never converge since can't be negative
    end
    
    TMeigs=zeros(max_modes,1); TEeigs=zeros(max_modes,1); hasConvergedTM=zeros(max_modes,2); hasConvergedTE=zeros(max_modes,2);
    N=0; TM_xaxis=[]; TE_xaxis=[];
    while 1
        if N<=30
            if N>=13 %here the TM converged in empty
                break
            end
            N=N+1;
        elseif N>30 & N <= 50
            N=N+3; disp(['N = ' num2str(N)])
        elseif N>50 & N <= 70
            N=N+10; disp(['N = ' num2str(N)])
        elseif N>70 & N <= 100
            if N>=72 %here the TE converged in empty
                break
            end
            N=N+15; disp(['N = ' num2str(N)])
        else
            disp('Necessary N for convergence exceeds 100. Aborting procedure.')
            break
        end        
        
        [TM,TE]=ComputeEigs(N, scenario);
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
    %convert to Number of Unknowns by #nodes/edges-#BCnodes/edges
    NoU_TE=@(N) N.^2-4*N+4;
    NoU_TM=@(N) N.*(N-1)+(2*N-1).*(N-1)-4*(N-1);
    
    %plot
    figure('Name','TM resonant frequencies k^2')
    hold on  
    plot(NoU_TM(TM_xaxis),TMeigs(:,TM_xaxis))
    %scatter(NoU_TM(hasConvergedTM(:,1)),hasConvergedTM(:,2),'filled')    
    %title('Error convergence of resonant frequencies as function of number of unknowns, \n for a cylinder completely filled with homogeneous (air) dielectric in the TM case','interpreter','latex')
    legend({'M_{edge}','G_{edge}'},'FontSize',8,'Location','best')
    xlabel('Number of unknowns','interpreter','latex')
    ylabel('Memory consumed in bytes','interpreter','latex')
    hold off
    
    figure('Name','TE resonant frequencies k^2')    
    hold on
    plot(NoU_TE(TE_xaxis),TEeigs(:,TE_xaxis))
    %scatter(NoU_TE(hasConvergedTE(:,1)),hasConvergedTE(:,2),'filled')        
    %title('Error convergence of resonant frequencies as function of number of unknowns, \n for a cylinder completely filled with homogeneous (air) dielectric in the TE case','interpreter','latex')
    legend({'M_{node}','G_{node}'},'FontSize',8,'Location','best')
    xlabel('Number of unknowns','interpreter','latex')
    ylabel('Memory consumed in bytes','interpreter','latex')    
    hold off
    
    figure('Name','TE resonant frequencies k^2')    
    hold on
    plot(NoU_TE(TE_xaxis),[TEeigs(:,TE_xaxis);TMeigs(:,TE_xaxis)])
    %scatter(NoU_TE(hasConvergedTE(:,1)),hasConvergedTE(:,2),'filled')        
    %title('Error convergence of resonant frequencies as function of number of unknowns, \n for a cylinder completely filled with homogeneous (air) dielectric in the TE case','interpreter','latex')
    legend({'M_{node}','G_{node}','M_{edge}','G_{edge}'},'FontSize',8,'Location','best')
    xlabel('Number of unknowns','interpreter','latex')
    ylabel('Memory consumed in bytes','interpreter','latex')    
    hold off
end
