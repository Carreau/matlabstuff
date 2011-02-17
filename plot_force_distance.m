fs=15

open('no_vasp_no_capping.mat');
plot(save_data(1).d,1e12*squeeze(save_data(1).f(1,1,:)))
xlabel('Bead distance [µm]','FontSize',fs+5);
ylabel('Force [pN]','FontSize',fs+5)
      
        set(gca,'FontSize',fs+5)
        set(gca,'LineWidth',2)
