plot(1:50,[mod(0.3*(1:50),1);mod(0.7*(1:50),1);mod(0.3*(1:50),1).*mod(0.7*(1:50),1)])
xlabel('$N$','interpreter','latex')
legend({'mod(0.3*N)','mod(0.7*N)','mod(0.3*N) \cdot mod(0.7*N)'},'FontSize',8,'Location','best')