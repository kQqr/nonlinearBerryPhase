tic

tolAlpha1 = 1e-7;
tolAlpha2 = 1e-5;
tolAlpha3 = 1e-3;
tolAlpha4 = 3e-3;

xAxisMax = 1.123;

for n80 = 1:2
    deltaX = 0.00125;
    if n80 == 1
        stdCutNum = 100;
    elseif n80 == 2
        stdCutNum = 1e10;
    end
    
    
    filename = ['eigM2Collect1_4_',num2str(1)];
    % saveas(gca,filename);
    load(filename);
    
    eigM2CollectTemp = zeros(2,size(eigM2Collect,2));
    label = 0;
    
    for eigM2List = [1 2 3 4 20 21 30 31 32 33 34 40 41 42]
        filename = ['eigM2Collect1_4_',num2str(eigM2List)];
        load(filename);
        
        for n = 1:size(eigM2Collect,1)
            if ~isnan(eigM2Collect(n,1))
                label = label+1;
                eigM2CollectTemp(label,1:size(eigM2Collect,2)) = eigM2Collect(n,1:size(eigM2Collect,2));
            end
        end
    end
    
    ampRec = zeros(2,1);
    labelAmp = 0;
    for n = 1:size(eigM2CollectTemp,1)
        flag = 1;
        for n1 = 1:labelAmp
            if eigM2CollectTemp(n,1) == ampRec(n1)
                flag = 0;
                break
            end
        end
        if flag == 1
            labelAmp = labelAmp+1;
            ampRec(labelAmp) = eigM2CollectTemp(n,1);
        end
    end
    ampRec = sort(ampRec);
    
    
    
    eigM2CollectTemp2 = zeros(2,size(eigM2Collect,2));
    label = 0;
    
    for n = 1:size(ampRec,1)
        amp = ampRec(n);
        for n1 = 1:size(eigM2CollectTemp,1)
            if eigM2CollectTemp(n1,1) == amp
                label = label+1;
                eigM2CollectTemp2(label,1:end) = eigM2CollectTemp(n1,1:end);
            end
        end
    end
    
    
    eigM2CollectMean_1_4 = zeros(1,size(eigM2CollectTemp2,2));
    for n = 1:size(ampRec,1)
        eigM2CollectAmp = zeros(1,size(eigM2CollectTemp2,2));
        labelTemp = 0;
        amp = ampRec(n);
        for n1 = 1:size(eigM2CollectTemp2,1)
            if eigM2CollectTemp2(n1,1) == amp
                if labelTemp <= stdCutNum
                    labelTemp = labelTemp+1;
                    eigM2CollectAmp(labelTemp,1:end) = eigM2CollectTemp2(n1,1:end);
                end
            end
        end
        eigM2CollectMean_1_4(n,1:end) = eigM2CollectAmp(1,1:end);
        eigM2CollectMean_1_4(n,14) = mean(eigM2CollectAmp(1:end,14));
        eigM2CollectMean_1_4(n,15) = std(eigM2CollectAmp(1:end,14));
        
        Pnum1 = 0;
        Pnum2 = 0;
        Pnum3 = 0;
        Pnum4 = 0;
        for n19 = 1:size(eigM2CollectAmp,1)
            if eigM2CollectAmp(n19,14) > tolAlpha1
                Pnum1 = Pnum1+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha2
                Pnum2 = Pnum2+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha3
                Pnum3 = Pnum3+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha4
                Pnum4 = Pnum4+1;
            end
        end
        eigM2CollectMean_1_4(n,16) = Pnum1/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_4(n,17) = Pnum2/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_4(n,18) = Pnum3/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_4(n,19) = Pnum4/size(eigM2CollectAmp,1);
    end
    
    % figure(111); clf;
    % hold on;
    % plot(eigM2CollectMean_1_4(:,1),eigM2CollectMean_1_4(:,14),'k.-')
    % for n = 1:size(eigM2CollectMean_1_4,1)
    %     plot([eigM2CollectMean_1_4(n,1) eigM2CollectMean_1_4(n,1)],[eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2],'k-')
    %     plot([eigM2CollectMean_1_4(n,1)-deltaX eigM2CollectMean_1_4(n,1)+deltaX],[eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2],'k-')
    %     plot([eigM2CollectMean_1_4(n,1)-deltaX eigM2CollectMean_1_4(n,1)+deltaX],[eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2],'k-')
    % end
    % %set(gca, 'YScale', 'log')
    % %axis([0 1.5 1e-8 5e-2])
    % hold off;
    
    eigM2Collect = eigM2CollectTemp;
    save('eigM2Collect1_4_0','eigM2Collect');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    eigM2CollectTemp = zeros(2,size(eigM2Collect,2));
    label = 0;
    
    for eigM2List = [1 3 4 5 6 7 8 9 10 11 12 13 15 16 30]
        filename = ['eigM2Collect1_5_',num2str(eigM2List)];
        load(filename);
        
        for n = 1:size(eigM2Collect,1)
            if ~isnan(eigM2Collect(n,1))
                label = label+1;
                eigM2CollectTemp(label,1:size(eigM2Collect,2)) = eigM2Collect(n,1:size(eigM2Collect,2));
            end
        end
    end
    
    ampRec = zeros(2,1);
    labelAmp = 0;
    for n = 1:size(eigM2CollectTemp,1)
        flag = 1;
        for n1 = 1:labelAmp
            if eigM2CollectTemp(n,1) == ampRec(n1)
                flag = 0;
                break
            end
        end
        if flag == 1
            labelAmp = labelAmp+1;
            ampRec(labelAmp) = eigM2CollectTemp(n,1);
        end
    end
    ampRec = sort(ampRec);
    
    
    
    eigM2CollectTemp2 = zeros(2,size(eigM2Collect,2));
    label = 0;
    
    for n = 1:size(ampRec,1)
        amp = ampRec(n);
        for n1 = 1:size(eigM2CollectTemp,1)
            if eigM2CollectTemp(n1,1) == amp
                label = label+1;
                eigM2CollectTemp2(label,1:end) = eigM2CollectTemp(n1,1:end);
            end
        end
    end
    
    
    eigM2CollectMean_1_5 = zeros(1,size(eigM2CollectTemp2,2));
    for n = 1:size(ampRec,1)
        eigM2CollectAmp = zeros(1,size(eigM2CollectTemp2,2));
        labelTemp = 0;
        amp = ampRec(n);
        for n1 = 1:size(eigM2CollectTemp2,1)
            if eigM2CollectTemp2(n1,1) == amp
                if labelTemp <= stdCutNum
                    labelTemp = labelTemp+1;
                    eigM2CollectAmp(labelTemp,1:end) = eigM2CollectTemp2(n1,1:end);
                end
            end
        end
        eigM2CollectMean_1_5(n,1:end) = eigM2CollectAmp(1,1:end);
        eigM2CollectMean_1_5(n,14) = mean(eigM2CollectAmp(1:end,14));
        eigM2CollectMean_1_5(n,15) = std(eigM2CollectAmp(1:end,14));
        
        Pnum1 = 0;
        Pnum2 = 0;
        Pnum3 = 0;
        Pnum4 = 0;
        for n19 = 1:size(eigM2CollectAmp,1)
            if eigM2CollectAmp(n19,14) > tolAlpha1
                Pnum1 = Pnum1+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha2
                Pnum2 = Pnum2+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha3
                Pnum3 = Pnum3+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha4
                Pnum4 = Pnum4+1;
            end
        end
        eigM2CollectMean_1_5(n,16) = Pnum1/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_5(n,17) = Pnum2/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_5(n,18) = Pnum3/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_5(n,19) = Pnum4/size(eigM2CollectAmp,1);

    end
    
    eigM2Collect = eigM2CollectTemp;
    save('eigM2Collect1_5_0','eigM2Collect');
    
%     figure(111); clf;
%     hold on;
%     plot(eigM2CollectMean_1_5(:,1),eigM2CollectMean_1_5(:,14),'k.-')
%     for n = 1:size(eigM2CollectMean_1_5,1)
%         plot([eigM2CollectMean_1_5(n,1) eigM2CollectMean_1_5(n,1)],[eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2],'k-')
%         plot([eigM2CollectMean_1_5(n,1)-deltaX eigM2CollectMean_1_5(n,1)+deltaX],[eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2],'k-')
%         plot([eigM2CollectMean_1_5(n,1)-deltaX eigM2CollectMean_1_5(n,1)+deltaX],[eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2],'k-')
%     end
%     %set(gca, 'YScale', 'log')
%     %axis([0 1.5 1e-8 5e-2])
%     hold off;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    eigM2CollectTemp = zeros(2,size(eigM2Collect,2));
    label = 0;
    
    for eigM2List = [1 2 3 4 5 6 7 8 9 10 11 12 13 20 21 30 31 32 33 34 35]
        filename = ['eigM2Collect1_7_',num2str(eigM2List)];
        load(filename);
        
        for n = 1:size(eigM2Collect,1)
            if ~isnan(eigM2Collect(n,1))
                label = label+1;
                eigM2CollectTemp(label,1:size(eigM2Collect,2)) = eigM2Collect(n,1:size(eigM2Collect,2));
            end
        end
    end
    
    ampRec = zeros(2,1);
    labelAmp = 0;
    for n = 1:size(eigM2CollectTemp,1)
        flag = 1;
        for n1 = 1:labelAmp
            if eigM2CollectTemp(n,1) == ampRec(n1)
                flag = 0;
                break
            end
        end
        if flag == 1
            labelAmp = labelAmp+1;
            ampRec(labelAmp) = eigM2CollectTemp(n,1);
        end
    end
    ampRec = sort(ampRec);
    
    
    
    eigM2CollectTemp2 = zeros(2,size(eigM2Collect,2));
    label = 0;
    
    for n = 1:size(ampRec,1)
        amp = ampRec(n);
        for n1 = 1:size(eigM2CollectTemp,1)
            if eigM2CollectTemp(n1,1) == amp
                label = label+1;
                eigM2CollectTemp2(label,1:end) = eigM2CollectTemp(n1,1:end);
            end
        end
    end
    
    
    eigM2CollectMean_1_7 = zeros(1,size(eigM2CollectTemp2,2));
    for n = 1:size(ampRec,1)
        eigM2CollectAmp = zeros(1,size(eigM2CollectTemp2,2));
        labelTemp = 0;
        amp = ampRec(n);
        for n1 = 1:size(eigM2CollectTemp2,1)
            if eigM2CollectTemp2(n1,1) == amp
                if labelTemp <= stdCutNum
                    labelTemp = labelTemp+1;
                    eigM2CollectAmp(labelTemp,1:end) = eigM2CollectTemp2(n1,1:end);
                end
            end
        end
        eigM2CollectMean_1_7(n,1:end) = eigM2CollectAmp(1,1:end);
        eigM2CollectMean_1_7(n,14) = mean(eigM2CollectAmp(1:end,14));
        eigM2CollectMean_1_7(n,15) = std(eigM2CollectAmp(1:end,14));
        
        Pnum1 = 0;
        Pnum2 = 0;
        Pnum3 = 0;
        Pnum4 = 0;
        for n19 = 1:size(eigM2CollectAmp,1)
            if eigM2CollectAmp(n19,14) > tolAlpha1
                Pnum1 = Pnum1+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha2
                Pnum2 = Pnum2+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha3
                Pnum3 = Pnum3+1;
            end
            if eigM2CollectAmp(n19,14) > tolAlpha4
                Pnum4 = Pnum4+1;
            end
        end
        eigM2CollectMean_1_7(n,16) = Pnum1/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_7(n,17) = Pnum2/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_7(n,18) = Pnum3/size(eigM2CollectAmp,1);
        eigM2CollectMean_1_7(n,19) = Pnum4/size(eigM2CollectAmp,1);

    end
    
    eigM2Collect = eigM2CollectTemp;
    save('eigM2Collect1_7_0','eigM2Collect');
    
    figure(111+n80); clf;
    hold on;
    
    H1 = plot([110 111],[110 111],'k.-','linewidth',1.5,'markersize',15);
    H2 = plot([110 111],[110 111],'b.-','linewidth',1.5,'markersize',15);
    H3 = plot([110 111],[110 111],'r.-','linewidth',1.5,'markersize',15);
    
    plot(eigM2CollectMean_1_4(:,1),eigM2CollectMean_1_4(:,14),'k.-')
    for n = 1:size(eigM2CollectMean_1_4,1)
        plot([eigM2CollectMean_1_4(n,1) eigM2CollectMean_1_4(n,1)],[eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2],'k-')
        plot([eigM2CollectMean_1_4(n,1)-deltaX eigM2CollectMean_1_4(n,1)+deltaX],[eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2],'k-')
        plot([eigM2CollectMean_1_4(n,1)-deltaX eigM2CollectMean_1_4(n,1)+deltaX],[eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2],'k-')
    end
    
    plot(eigM2CollectMean_1_5(:,1),eigM2CollectMean_1_5(:,14),'b.-')
    for n = 1:size(eigM2CollectMean_1_5,1)
        plot([eigM2CollectMean_1_5(n,1) eigM2CollectMean_1_5(n,1)],[eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2],'b-')
        plot([eigM2CollectMean_1_5(n,1)-deltaX eigM2CollectMean_1_5(n,1)+deltaX],[eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2],'b-')
        plot([eigM2CollectMean_1_5(n,1)-deltaX eigM2CollectMean_1_5(n,1)+deltaX],[eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2],'b-')
    end
    
    plot(eigM2CollectMean_1_7(:,1),eigM2CollectMean_1_7(:,14),'r.-')
    for n = 1:size(eigM2CollectMean_1_7,1)
        plot([eigM2CollectMean_1_7(n,1) eigM2CollectMean_1_7(n,1)],[eigM2CollectMean_1_7(n,14)-eigM2CollectMean_1_7(n,15)/2 eigM2CollectMean_1_7(n,14)+eigM2CollectMean_1_7(n,15)/2],'r-')
        plot([eigM2CollectMean_1_7(n,1)-deltaX eigM2CollectMean_1_7(n,1)+deltaX],[eigM2CollectMean_1_7(n,14)-eigM2CollectMean_1_7(n,15)/2 eigM2CollectMean_1_7(n,14)-eigM2CollectMean_1_7(n,15)/2],'r-')
        plot([eigM2CollectMean_1_7(n,1)-deltaX eigM2CollectMean_1_7(n,1)+deltaX],[eigM2CollectMean_1_7(n,14)+eigM2CollectMean_1_7(n,15)/2 eigM2CollectMean_1_7(n,14)+eigM2CollectMean_1_7(n,15)/2],'r-')
    end
    
    
    %set(gca, 'YScale', 'log')
    axis([0.7 xAxisMax 0 0.025])
    hold off;
    legend([H1 H2 H3],{'$q=\pi/2$','$q=2\pi/5$','$q=2\pi/7$'},'Location','northwest','FontSize',19,'Interpreter','latex');
    xticks([0.7 0.9 1.1])
    xlabel('$A$','Interpreter','latex')
    %yticks([1e-6 1e-4 1e-2 1e0])
    %yticklabels({'-6','-4','-2','0'})
    yticks([0.01 0.02])
    ylabel('$\langle \max\{{\rm Re\,}\alpha_n\}\rangle$','Interpreter','latex')
    % set(gca, 'YScale', 'log')
    set(gca,'linewidth',2)
    set(gca,'FontSize',19)
    
    
    
    
    figure(1110+n80); clf;
    hold on;
    
    H1 = plot([110 111],[110 111],'k.-','linewidth',1.5,'markersize',15);
    H2 = plot([110 111],[110 111],'b.-','linewidth',1.5,'markersize',15);
    H3 = plot([110 111],[110 111],'r.-','linewidth',1.5,'markersize',15);
    
    plot(eigM2CollectMean_1_4(:,1),eigM2CollectMean_1_4(:,14),'k.-')
    for n = 1:size(eigM2CollectMean_1_4,1)
        plot([eigM2CollectMean_1_4(n,1) eigM2CollectMean_1_4(n,1)],[eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2],'k-')
        plot([eigM2CollectMean_1_4(n,1)-deltaX eigM2CollectMean_1_4(n,1)+deltaX],[eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)-eigM2CollectMean_1_4(n,15)/2],'k-')
        plot([eigM2CollectMean_1_4(n,1)-deltaX eigM2CollectMean_1_4(n,1)+deltaX],[eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2 eigM2CollectMean_1_4(n,14)+eigM2CollectMean_1_4(n,15)/2],'k-')
    end
    
    plot(eigM2CollectMean_1_5(:,1),eigM2CollectMean_1_5(:,14),'b.-')
    for n = 1:size(eigM2CollectMean_1_5,1)
        plot([eigM2CollectMean_1_5(n,1) eigM2CollectMean_1_5(n,1)],[eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2],'b-')
        plot([eigM2CollectMean_1_5(n,1)-deltaX eigM2CollectMean_1_5(n,1)+deltaX],[eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)-eigM2CollectMean_1_5(n,15)/2],'b-')
        plot([eigM2CollectMean_1_5(n,1)-deltaX eigM2CollectMean_1_5(n,1)+deltaX],[eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2 eigM2CollectMean_1_5(n,14)+eigM2CollectMean_1_5(n,15)/2],'b-')
    end
    
    plot(eigM2CollectMean_1_7(:,1),eigM2CollectMean_1_7(:,14),'r.-')
    for n = 1:size(eigM2CollectMean_1_7,1)
        plot([eigM2CollectMean_1_7(n,1) eigM2CollectMean_1_7(n,1)],[eigM2CollectMean_1_7(n,14)-eigM2CollectMean_1_7(n,15)/2 eigM2CollectMean_1_7(n,14)+eigM2CollectMean_1_7(n,15)/2],'r-')
        plot([eigM2CollectMean_1_7(n,1)-deltaX eigM2CollectMean_1_7(n,1)+deltaX],[eigM2CollectMean_1_7(n,14)-eigM2CollectMean_1_7(n,15)/2 eigM2CollectMean_1_7(n,14)-eigM2CollectMean_1_7(n,15)/2],'r-')
        plot([eigM2CollectMean_1_7(n,1)-deltaX eigM2CollectMean_1_7(n,1)+deltaX],[eigM2CollectMean_1_7(n,14)+eigM2CollectMean_1_7(n,15)/2 eigM2CollectMean_1_7(n,14)+eigM2CollectMean_1_7(n,15)/2],'r-')
    end

    %set(gca, 'YScale', 'log')
    axis([0 xAxisMax 0 0.025])
    hold off;
    legend([H1 H2 H3],{'$q=\pi/2$','$q=2\pi/5$','$q=2\pi/7$'},'Location','northwest','FontSize',19,'Interpreter','latex');
    xticks([0 0.5 1])
    xlabel('$A$','Interpreter','latex')
    %yticks([1e-6 1e-4 1e-2 1e0])
    %yticklabels({'-6','-4','-2','0'})
    yticks([0.01 0.02])
    ylabel('$\langle \max\{{\rm Re\,}\alpha_n\}\rangle$','Interpreter','latex')
    % set(gca, 'YScale', 'log')
    set(gca,'linewidth',2)
    set(gca,'FontSize',19)
    
    
    
    
    
    figure(11100+n80); clf;
    hold on;
    
    H1 = plot([110 111],[110 111],'k.-','linewidth',1.5,'markersize',15);
    H2 = plot([110 111],[110 111],'b.-','linewidth',1.5,'markersize',15);
    H3 = plot([110 111],[110 111],'r.-','linewidth',1.5,'markersize',15);
    
    plot(eigM2CollectMean_1_4(:,1),eigM2CollectMean_1_4(:,18),'k.-')
    plot(eigM2CollectMean_1_5(:,1),eigM2CollectMean_1_5(:,18),'b.-')
    plot(eigM2CollectMean_1_7(:,1),eigM2CollectMean_1_7(:,18),'r.-')

    %set(gca, 'YScale', 'log')
    axis([0 xAxisMax 0 1])
    hold off;
    legend([H1 H2 H3],{'$q=\pi/2$','$q=2\pi/5$','$q=2\pi/7$'},'Location','northwest','FontSize',19,'Interpreter','latex');
    xticks([0 0.5 1])
    xlabel('$A$','Interpreter','latex')
    %yticks([1e-6 1e-4 1e-2 1e0])
    %yticklabels({'-6','-4','-2','0'})
    yticks([0 1])
    ylabel('$P(\langle \max\{{\rm Re\,}\alpha_n\}\rangle>10^{-3})$','Interpreter','latex')
    % set(gca, 'YScale', 'log')
    set(gca,'linewidth',2)
    set(gca,'FontSize',19)
    
end





toc










