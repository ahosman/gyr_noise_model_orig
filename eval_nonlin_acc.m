function eval_nonlin_acc(obj)
    set(0,...
        'defaultFigureColor','w',...
        'DefaultFigureWindowStyle','docked', ...
        'defaultTextInterpreter','tex', ...
        'defaultAxesFontSize',14, ...
        'defaultLegendInterpreter','tex', ...
        'defaultLineLineWidth', 1);
    format long g
    h1 = figure;
    h2 = figure;
    h3 = figure;
    h4 = figure;
    
    handles_fig2 = [];
    handles_fig3 = [];
    handles_fig4 = [];
        
    figure(h1);
    g = -64:0.1:64;
    C = obj.acc_mems(obj.Sg, g, obj.beta, 1);
    
    %% Original settings
    fprintf("\n#######################################");
    fprintf("\n        Original Settings              ");
    fprintf("\n#######################################");
    color = 'r';
    gtrim = 32;
    figure(h2);
    [hp1, Vcv_32g, Vcv_8g, acc_afe_gain_corr_32g, acc_afe_gain_corr_8g, acc_cic_gain_corr, acc_trm_dgain_s0] = ...
        obj.c2v_acc_afe(obj.Sg, g, C, obj.vsat, gtrim, 'original', obj.vref, obj.vcv_max, 1, color);
    handles_fig2 = [handles_fig2, hp1];

    figure(h3);
    hd11 = obj.acc_dp('32g', g, obj.vref, 1, color, Vcv_32g, acc_afe_gain_corr_32g, acc_cic_gain_corr, acc_trm_dgain_s0, 'original');
    handles_fig3 = [handles_fig3, hd11];

    figure(h4);
    hd12 = obj.acc_dp('8g' , g, obj.vref, 1, color, Vcv_8g , acc_afe_gain_corr_8g , acc_cic_gain_corr, acc_trm_dgain_s0);
    handles_fig4 = [handles_fig4, hd12];
       
    %% Option 1(b): Double-range extension 0.85V @ 40g
    fprintf("\n#######################################");
    fprintf("\n        Double-Range @ 40g             ");
    fprintf("\n#######################################");
    color = 'b';
    gtrim = 40;
    
    figure(h2);
    [hp2, Vcv_32g, Vcv_8g, acc_afe_gain_corr_32g, acc_afe_gain_corr_8g, acc_cic_gain_corr, acc_trm_dgain_s0] = ...
        obj.c2v_acc_afe(obj.Sg, g, C, obj.vsat, gtrim, 'original', obj.vref, obj.vcv_max, 1, color);
    handles_fig2 = [handles_fig2, hp2];
    
    figure(h3);
    hd31 = obj.acc_dp('32g', g, obj.vref, 1, color, Vcv_32g, acc_afe_gain_corr_32g, acc_cic_gain_corr, acc_trm_dgain_s0, 'dbl_rng_40g',0);
    handles_fig3 = [handles_fig3, hd31];
    
    figure(h4);
    hd32 = obj.acc_dp('8g',  g, obj.vref, 1, color, Vcv_8g , acc_afe_gain_corr_8g , acc_cic_gain_corr, acc_trm_dgain_s0);
    handles_fig4 = [handles_fig4, hd32];
        
    figure(h2);
    xline(+32, 'r--', 'LineWidth', 1);
    xline(-32, 'r--', 'LineWidth', 1);
    xline( +8, 'r--', 'LineWidth', 1);
    xline( -8, 'r--', 'LineWidth', 1);
    yline(+0.85, 'g--', 'LineWidth',1);
    yline(-0.85, 'g--', 'LineWidth',1);
    grid on;
    xlabel('Input acceleration (g)');
    ylabel('C2V Output Voltage (V)');
    title(sprintf('C2V Output Voltage (vsat = %.2f)',obj.vsat));
    ylim([-1.5 1.5])
    legend(handles_fig2(1,:), {'original','double-range:0.85@40g'}, ...
        'Position',[0.5546875 0.381027667984189 0.156305625000005 0.0988142292490124]);
    
    figure(h3);
    xline(+32, 'r--', 'LineWidth', 1);
    xline(-32, 'r--', 'LineWidth', 1);
    xlabel('Input acceleration (g)');
    ylabel('Digital Output acceleration (g)');
    title(sprintf('Digital Output-32g'));
    xlim([-40 40])
    grid on;
    legend(handles_fig3(1,:), {'original','double-range:0.85@40g'}, ...
        'Position',[0.562506604347373 0.279533571808121 0.148708730468748 0.135948687747036]);
    
    figure(h4);
    xline(+8, 'r--', 'LineWidth', 1);
    xline(-8, 'r--', 'LineWidth', 1);
    xlabel('Input acceleration (g)');
    ylabel('Digital Output acceleration (g)');
    title(sprintf('Digital Output-8g'));
    xlim([-10 10])
    grid on;
    legend(handles_fig4(1,:), {'original','double-range:0.85@40g'}, ...
        'Position',[0.562506604347373 0.279533571808121 0.148708730468748 0.135948687747036]);
    fprintf("\n");
end
