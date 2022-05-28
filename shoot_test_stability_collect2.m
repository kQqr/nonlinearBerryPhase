tic

rng("shuffle");

load('eigM2Collect1_7_0');
eigM2CollectTemp = eigM2Collect;

dAcut = 0.025;%0.0000001;
fftNum = 100; 
randLabelMax = 500;
startFFTloop = 1;
randA = 1e-3;

tol1 = 1e-6;
wRecSwitch = 0;
figSee = 0;

% basic physical parameters of the 1D chain
load('NNM1_7pp');
%load('NNM2_5pp');
NNMtemp = NNM;

ASortTemp = nan(size(NNM,1)/4,1);
for n1 = 1:size(NNM,1)/4
    ASortTemp(n1) = NNM(4*(n1-1)+2,1);
end
[ASort,ASortOrder] = sort(ASortTemp);
for n1 = 1:size(NNM,1)/4
    NNM(4*(n1-1)+1,:) = NNMtemp(4*(ASortOrder(n1)-1)+1,:);
    NNM(4*(n1-1)+2,:) = NNMtemp(4*(ASortOrder(n1)-1)+2,:);
    NNM(4*(n1-1)+3,:) = NNMtemp(4*(ASortOrder(n1)-1)+3,:);
    NNM(4*(n1-1)+4,:) = NNMtemp(4*(ASortOrder(n1)-1)+4,:);
end

wRec = nan(1000,1); % record the omega during the process of calculating the NNMs
wPer = nan(1000,1);
wRecUp = size(NNM,1)/4;

% for n12 = 1:wRecUp
%     wRec(n12) = NNM(size(NNM,1)-4-4*(n12-1)+2,2);
%     if n12 >= 2
%         wPer(n12-1) = (wRec(n12)-wRec(n12-1))/wRec(n12-1);
%     end
% end
for n12 = 1:wRecUp
    wRec(n12) = NNM(4*(n12-1)+2,2);
    if n12 >= 2
        wPer(n12-1) = (wRec(n12)-wRec(n12-1))/wRec(n12-1);
    end
end
% NNMlabel = 141; %size(NNM,1)/4;
% NNMlabel = floor(size(NNM,1)/4);

eigM2Collect = nan(randLabelMax*size(NNM,1)/4,20);
NNMlabelRand = 0;

Apre = 0;
for NNMlabel = 169:206
    Anow = NNM(4*(NNMlabel-1)+2,1);
    flag0 = 0;
    for n78 = 1:size(eigM2CollectTemp,1)
        if abs(eigM2CollectTemp(n78,1)-Anow) < 1e-7
            flag0 = 1;
            break
        end
    end
    if Anow-Apre < dAcut
        continue
    elseif flag0 == 1
        continue
    else
        Apre = NNM(4*(NNMlabel-1)+2,1);
        for randLabel = 1:randLabelMax
            
            NNMlabelRand = NNMlabelRand+1;
            
            nStart = size(NNM,1);%size(NNM,1);
            
            e0 = 1.5;%2;
            c1 = 0.25;%0.1;%1
            c2 = NNM(2,5)*0.37;%0.35;%2
            
            d = 0;
            d1 = 0.22;%0.22;
            d2 = NNM(2,5)*0.02;%0.02;%0.02;
            
            f1 = d1;
            f2 = d2;
            
            Ac = sqrt(-4*(abs(c2)-c1)/3/(abs(d2)-d1));
            
            ucNum = size(NNM,2)/4;
            sw = NNM(4*(1-1)+2,4);
            scd = NNM(4*(1-1)+2,5);
            tol2 = 1e-2;
            tol3 = 1e-4;
            tol4 = 1e-4;
            
            q = NNM(4*(1-1)+2,7);
            DOF = 4*ucNum;
            timeN = 1000;
            critLoopNum = 10000;
            nu = 1/4;
            
            % initialize the NNM using shooting method
            zk = zeros(timeN+1,DOF);
            gidt = zeros(7,DOF);
            CorrectionList = zeros(DOF,2);
            CorrectionLabel = 0;
            
            zfft1 = nan(timeN,1);
            
            maxzk = Ac*ones(DOF,1); % the peak of the disp
            shiftzk = zeros(timeN+1,DOF); % shift the displacement such that the time starts from the maximum position
            shiftLabel = zeros(DOF,1); % record the shifting label
            shiftT = zeros(DOF,1); % record the shifting time
            IntraLabel = zeros(DOF,1); % record the intracell label shift
            IntraPhase = zeros(DOF,1); % record the intracell phase shift
            InterLabel = zeros(DOF,1); % record the intercell label shift
            InterPhase = zeros(DOF,1); % record the intercell phase shift
            
            zk(1,1:DOF) = NNM(4*(NNMlabel-1)+1,1:DOF);
            for n = 1:DOF
                zk(1,n) = zk(1,n)+randA*rand(1);
            end
            Tk = 2*pi/NNM(4*(NNMlabel-1)+2,2);
            dTk = Tk/timeN;
            
            loopA = wRecUp;
            loopA1 = wRecUp;
            
            errorRec = nan(critLoopNum+1,1);
            
            M3 = zeros(DOF,DOF);
            
            eigM2 = zeros(2,1);
            
            labelz = 0;
            for loopNum = 1:fftNum
                if loopNum ~= 1
                    zk(1,1:DOF) = zk(timeN+1,1:DOF);
                end
                
                M2 = zeros(DOF,DOF);
                
                for n = 1:timeN
                    
                    Midt = zeros(DOF,DOF,7);
                    yi = zeros(7,DOF);
                    
                    for n20 = 1:7
                        if n20 == 1
                            yi(1,:) = zk(n,:);
                        elseif n20 == 2
                            yi(2,:) = yi(1,:)+gidt(1,:)*nu;
                        elseif n20 == 3
                            yi(3,:) = yi(1,:)+((4*nu-1)*gidt(1,:)+gidt(2,:))/8/nu;
                        elseif n20 == 4
                            yi(4,:) = yi(1,:)+((10*nu-2)*gidt(1,:)+2*gidt(2,:)+8*nu*gidt(3,:))/27/nu;
                        elseif n20 == 5
                            yi(5,:) = yi(1,:)+(-(77*nu-56+(17*nu-8)*sqrt(21))*gidt(1,:)-8*(7+sqrt(21))*gidt(2,:)+48*(7+sqrt(21))*nu*gidt(3,:)-3*(21+sqrt(21))*nu*gidt(4,:))/nu/392;
                        elseif n20 == 6
                            yi(6,:) = yi(1,:)+(-5*(287*nu-56-(59*nu-8)*sqrt(21))*gidt(1,:)-40*(7-sqrt(21))*gidt(2,:)+320*sqrt(21)*nu*gidt(3,:)+3*(21-121*sqrt(21))*nu*gidt(4,:)+392*(6-sqrt(21))*nu*gidt(5,:))/nu/1960;
                        elseif n20 == 7
                            yi(7,:) = yi(1,:)+(15*(30*nu-8-7*nu*sqrt(21))*gidt(1,:)+120*gidt(2,:)-40*(5+7*sqrt(21))*nu*gidt(3,:)+63*(2+3*sqrt(21))*nu*gidt(4,:)-14*(49-9*sqrt(21))*nu*gidt(5,:)+70*(7+sqrt(21))*nu*gidt(6,:))/nu/180;
                        end
                        
                        for n1 = 1:ucNum
                            if n1 == 1
                                n1temp = n1+ucNum;
                            else
                                n1temp = n1;
                            end
                            if n1 == ucNum
                                n2temp = n1-ucNum;
                            else
                                n2temp = n1;
                            end
                            gidt(n20,4*n1-3) = +dTk*(+e0*yi(n20,4*n1-2)+d*yi(n20,4*n1-2)^3+c1*yi(n20,4*n1+0)+f1*yi(n20,4*n1+0)^3+c2*yi(n20,4*n1temp-4)+f2*yi(n20,4*n1temp-4)^3);
                            gidt(n20,4*n1-2) = -dTk*(+e0*yi(n20,4*n1-3)+d*yi(n20,4*n1-3)^3+c1*yi(n20,4*n1-1)+d1*yi(n20,4*n1-1)^3+c2*yi(n20,4*n1temp-5)+d2*yi(n20,4*n1temp-5)^3);
                            gidt(n20,4*n1-1) = +dTk*(+e0*yi(n20,4*n1-0)+d*yi(n20,4*n1-0)^3+c1*yi(n20,4*n1-2)+f1*yi(n20,4*n1-2)^3+c2*yi(n20,4*n2temp+2)+f2*yi(n20,4*n2temp+2)^3);
                            gidt(n20,4*n1-0) = -dTk*(+e0*yi(n20,4*n1-1)+d*yi(n20,4*n1-1)^3+c1*yi(n20,4*n1-3)+d1*yi(n20,4*n1-3)^3+c2*yi(n20,4*n2temp+1)+d2*yi(n20,4*n2temp+1)^3);
                            
                            Midt(4*n1-3,4*n1temp-4,n20) = dTk*(+c2+3*f2*yi(n20,4*n1temp-4)^2);
                            Midt(4*n1-3,4*n1-2,n20)     = dTk*(+e0+3*d*yi(n20,4*n1-2)^2);
                            Midt(4*n1-3,4*n1-0,n20)     = dTk*(+c1+3*f1*yi(n20,4*n1-0)^2);
                            
                            Midt(4*n1-2,4*n1temp-5,n20) = dTk*(-c2-3*d2*yi(n20,4*n1temp-5)^2);
                            Midt(4*n1-2,4*n1-3,n20)     = dTk*(-e0-3*d*yi(n20,4*n1-3)^2);
                            Midt(4*n1-2,4*n1-1,n20)     = dTk*(-c1-3*d1*yi(n20,4*n1-1)^2);
                            
                            Midt(4*n1-1,4*n1-2,n20)     = dTk*(+c1+3*f1*yi(n20,4*n1-2)^2);
                            Midt(4*n1-1,4*n1-0,n20)     = dTk*(+e0+3*d*yi(n20,4*n1-0)^2);
                            Midt(4*n1-1,4*n2temp+2,n20) = dTk*(+c2+3*f2*yi(n20,4*n2temp+2)^2);
                            
                            Midt(4*n1-0,4*n1-3,n20)     = dTk*(-c1-3*d1*yi(n20,4*n1-3)^2);
                            Midt(4*n1-0,4*n1-1,n20)     = dTk*(-e0-3*d*yi(n20,4*n1-1)^2);
                            Midt(4*n1-0,4*n2temp+1,n20) = dTk*(-c2-3*d2*yi(n20,4*n2temp+1)^2);
                            
                        end
                    end
                    
                    zk(n+1,:) = zk(n,:)+(9*gidt(1,:)+64*gidt(3,:)+49*gidt(5,:)+49*gidt(6,:)+9*gidt(7,:))/180;
                    M2 = M2+(9*Midt(:,:,1)+64*Midt(:,:,3)+49*Midt(:,:,5)+49*Midt(:,:,6)+9*Midt(:,:,7))/180;
                end
                if loopNum >= startFFTloop
                    labelz = labelz+1;
                    zfft1(timeN*(labelz-1)+1:timeN*labelz) = zk(1:timeN,1);
                end
                
                M3 = M3+M2;
                
                outM = expm(M2)-eye(DOF,DOF);
                
                % eigM2(loopNum) = max(abs(log(abs(eig(expm(M2))))));
                eigM2(loopNum) = max(log(abs(eig(expm(M2)))));
                
                %     if loopNum == 1 && figSee == 1
                %         [maxzk1,maxzk1Label] = max(zk(1:timeN+1,4*(1-1)+1));
                %         cosCompare = zeros(timeN+1,1);
                %         for n2 = 1:timeN+1
                %             cosCompare(n2) = maxzk1*cos((n2-maxzk1Label)/timeN*2*pi);
                %         end
                %         figure(2); clf;
                %         maxzk1 = max(zk(:,1));
                %         hold on;
                %         H1 = plot(dTk*(1:timeN),zk(1:timeN,4*(1-1)+1),'r-','linewidth',3);
                %         H2 = plot(dTk*(1:timeN),cosCompare(1:timeN),'k-.','linewidth',3);
                %         H3 = plot(dTk*(1:timeN),zk(1:timeN,4*(1+1-1)+1),'b-','linewidth',3);
                %         axis([0 Tk -maxzk1 maxzk1])
                %         hold off;
                %         legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','$A\cos(\omega t)$','${\rm Re\,}\Psi_2^{(1)}(t)$'},'Location','east','FontSize',16,'Interpreter','latex');
                %         set(gca,'linewidth',2)
                %         set(gca,'FontSize',19)
                %         xlabel('$t$','Interpreter','latex')
                %         ylabel('${\rm Re\,}\Psi_n^{(1)}(t)$','Interpreter','latex')
                %         xticks([0 1.5 3])
                %         yticks([-1.5 0 1.5])
                %     end
                
            end
            
            % zfft1(timeN*(fftNum-1)+2)-zfft1(2)
            
            n1 = 1;
            [maxzk1,maxzk1Label] = max(zk(1:timeN+1,4*(n1-1)+1));
            cosCompare = zeros(timeN+1,1);
            for n2 = 1:timeN+1
                cosCompare(n2) = maxzk1*cos((n2-maxzk1Label)/timeN*2*pi);
            end
            if figSee == 1
                figure(1); clf;
                maxzk1 = max(zk(:,1));
                hold on;
                H1 = plot(dTk*(1:timeN),zk(1:timeN,4*(n1-1)+1),'r-','linewidth',3);
                H2 = plot(dTk*(1:timeN),cosCompare(1:timeN),'k-.','linewidth',3);
                H3 = plot(dTk*(1:timeN),zk(1:timeN,4*(n1+1-1)+3),'b-','linewidth',3);
                axis([0 Tk -maxzk1 maxzk1])
                hold off;
                legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','$A\cos(\omega t)$','${\rm Re\,}\Psi_2^{(2)}(t)$'},'Location','east','FontSize',16,'Interpreter','latex');
                set(gca,'linewidth',2)
                set(gca,'FontSize',19)
                xlabel('$t$','Interpreter','latex')
                ylabel('${\rm Re\,}\Psi_n^{(1)}(t)$','Interpreter','latex')
                xticks([0 1.5 3])
                yticks([-1.5 0 1.5])
                
                
                
                
                
                
                
                tol6 = 1e-3;
                P1x = 0;
                Q1x = 0;
                P2x = 0;
                Q2x = 0;
                
                FreqMeasure = 2*pi/dTk/size(zfft1,1)*(0:size(zfft1,1)); % the ruler to measure the spectrum: in unit of freqUnit
                
                P2x = 1/size(zfft1,1)*abs(fft(zfft1(1:end)));
                P1x = P2x;
                P1x = P2x(1:floor(size(P2x,1)/2));
                P1x(2:end-1) = 2*P1x(2:end-1);
                
                [maxP1x,maxP1xLabel] = max(P1x);
                fftLabelLimit = floor(maxP1xLabel*16);
                
                figure(9); clf;
                hold on;
                
                %     bsBot = e0+c1-c2;
                %     bsTop = e0-c1+c2;
                %     % X = [wNonlinear(1,1); WNonlinear(1,1); WNonlinear(1,1); wNonlinear(1,1)];
                %     X = [bsBot; bsTop; bsTop; bsBot];
                %     Y = [0; 0; 1; 1];
                %     h = fill(X,Y,gold1);
                %     set(h,'EdgeColor','none')
                %     alpha(alphaValue)
                
                H2 = plot(FreqMeasure(2:fftLabelLimit),P1x(2:fftLabelLimit),'b-','linewidth',3);
                % H1 = plot(FreqMeasure(1:fftLabelLimit),sqrt(2)*P1x(1:fftLabelLimit)/sqrt(2)/max(P1x(1:fftLabelLimit)),'b-','linewidth',LineWidth,'color',brown);
                % plot([+c1-c2 +c1-c2],[0 3e-4],'k--','linewidth',LineWidth)
                % plot([-c1+c2 -c1+c2],[0 3e-4],'k--','linewidth',LineWidth)
                hold off;
                set(gca, 'YScale', 'log')
                legend([H2],{'${\rm Re\,}\Psi_1^{(1)}(\omega)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
                set(gca,'linewidth',2)
                set(gca,'FontSize',19)
                
                xlabel('$\omega$','Interpreter','latex')
                ylabel('${\rm Re\,}\Psi_1^{(1)}(\omega)$','Interpreter','latex')
                xticks([0 10 20])
                yticks([1e-6 1e-3 1])
                axis([0 20 1e-6 maxzk1*1.1])
                
                %     filename = ['fft',num2str(fLabel),'.png'];
                %     saveas(gca,filename);
            end
            
            
            
            
            
            
            errorRec2 = errorRec;
            
            
            
            for n1 = 1:ucNum
                %     figure(n1); clf;
                %     hold on;
                %     plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+1))
                %     plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+2))
                %     plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+3))
                %     plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+4))
                %     axis([0 Tk -maxzk(1) maxzk(1)])
                %     hold off;
            end
            
            maxzk = zeros(DOF,1); % the peak of the disp
            shiftzk = zeros(timeN+1,DOF); % shift the displacement such that the time starts from the maximum position
            shiftLabel = zeros(DOF,1); % record the shifting label
            shiftT = zeros(DOF,1); % record the shifting time
            IntraLabel = zeros(DOF,1); % record the intracell label shift
            IntraPhase = zeros(DOF,1); % record the intracell phase shift
            InterLabel = zeros(DOF,1); % record the intercell label shift
            InterPhase = zeros(DOF,1); % record the intercell phase shift
            for n1 = 1:DOF
                [maxzk(n1),startLabel] = max(zk(1:timeN+1,n1));
                shiftLabel(n1) = startLabel;
                shiftT(n1) = dTk*shiftLabel(n1);
                nTemp = 0;
                for n = [startLabel:timeN+1 1:startLabel-1]
                    nTemp = nTemp+1;
                    shiftzk(nTemp,n1) = zk(n,n1);
                end
            end
            for n1 = 1:ucNum
                for n22 = 1:4
                    IntraLabel(4*(n1-1)+n22) = mod((shiftLabel(4*(n1-1)+n22)-shiftLabel(4*(n1-1)+1)),timeN);
                    %         if abs(IntraLabel(4*(n1-1)+n22)-timeN)/timeN < tol2 || abs(IntraLabel(4*(n1-1)+n22))/timeN < tol2
                    %             IntraLabel(4*(n1-1)+n22) = 0;
                    %         end
                    IntraPhase(4*(n1-1)+n22) = mod((shiftT(4*(n1-1)+n22)-shiftT(4*(n1-1)+1))/Tk,1);
                    %         if abs(IntraPhase(4*(n1-1)+n22)-1)/1 < tol3 || abs(IntraLabel(4*(n1-1)+n22))/1 < tol3
                    %             IntraPhase(4*(n1-1)+n22) = 0;
                    %         end
                    if n1 == ucNum
                        n2temp = n1-ucNum;
                    else
                        n2temp = n1;
                    end
                    InterLabel(4*(n1-1)+n22) = mod((shiftLabel(4*(n2temp-0)+n22)-shiftLabel(4*(n1-1)+n22)),timeN);
                    %         if abs(InterLabel(4*(n1-1)+n22)-timeN)/timeN < tol2 || abs(InterLabel(4*(n1-1)+n22))/timeN < tol2
                    %             InterLabel(4*(n1-1)+n22) = 0;
                    %         end
                    InterPhase(4*(n1-1)+n22) = mod((shiftT(4*(n2temp-0)+n22)-shiftT(4*(n1-1)+n22))/Tk,1);
                    %         if abs(InterPhase(4*(n1-1)+n22)-1)/1 < tol3 || abs(InterPhase(4*(n1-1)+n22))/1 < tol3
                    %             InterPhase(4*(n1-1)+n22) = 0;
                    %         end
                end
            end
            
            disp(IntraPhase(3))
            
            
            %         figure(3); clf;
            %         hold on;
            %
            %         H1 = plot(dTk*(1:size(zfft1,1))/Tk,zfft1(:),'b-','linewidth',0.5);
            %
            %         axis([0 dTk*size(zfft1,1)/Tk -1.6 1.6])
            %         hold off;
            %         set(gca,'linewidth',2)
            %         set(gca,'FontSize',19)
            %         xticks([0 100 200])
            %         yticks([-1.5 0 1.5])
            %         % ylabel('${\rm Re\,}\Psi_1^{(1)}(t)$','Interpreter','latex')
            %
            %         %legend([H1],{'${\rm Re\,}\Psi_1^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
            %         set(gca,'linewidth',2)
            %         set(gca,'FontSize',19)
            %         xlabel('$t/T$','Interpreter','latex')
            %
            %         % yticklabels({'-A_c','-0.5','0','0.5','A_c'})
            %         % axis([0 500 -Ac*1.05 Ac*1.05])
            
            
            % eigM3 = max(abs(log(abs(eig(expm(M3))))));
            
            % save('NNMtemp','NNM');
            
            % max(eigM2)
            % mean(eigM2)
            % std(eigM2)
            
            eigM2Collect(NNMlabelRand,1)  = NNM(4*(NNMlabel-1)+2,1); % A
            eigM2Collect(NNMlabelRand,2)  = NNM(4*(NNMlabel-1)+2,2); % w
            eigM2Collect(NNMlabelRand,3)  = NNM(4*(NNMlabel-1)+2,3); % error of shooting method
            eigM2Collect(NNMlabelRand,4)  = NNM(4*(NNMlabel-1)+2,4); % sw
            eigM2Collect(NNMlabelRand,5)  = NNM(4*(NNMlabel-1)+2,5); % scd
            eigM2Collect(NNMlabelRand,6)  = NNM(4*(NNMlabel-1)+2,6); % ucNum
            eigM2Collect(NNMlabelRand,7)  = NNM(4*(NNMlabel-1)+2,7); % q
            eigM2Collect(NNMlabelRand,8)  = NNM(4*(NNMlabel-1)+2,8); % c1
            eigM2Collect(NNMlabelRand,9)  = NNM(4*(NNMlabel-1)+2,9); % c2
            eigM2Collect(NNMlabelRand,10) = NNM(4*(NNMlabel-1)+2,10); % d1
            eigM2Collect(NNMlabelRand,11) = NNM(4*(NNMlabel-1)+2,11); % d2
            eigM2Collect(NNMlabelRand,12) = NNM(4*(NNMlabel-1)+2,12); % Ac
            eigM2Collect(NNMlabelRand,13) = fftNum; %
            eigM2Collect(NNMlabelRand,14) = mean(eigM2); %
            eigM2Collect(NNMlabelRand,15) = std(eigM2); %
            
            % min(errorRec2)
        end
    end
    
    save('eigM2Collect1_7_35','eigM2Collect');
    
end

% mean(eigM2Collect(1:randLabelMax,14))
% std(eigM2Collect(1:randLabelMax,14))

eigM2CollectMean = nan(size(eigM2Collect,1)/randLabelMax,size(eigM2Collect,2));
for n = 1:size(eigM2Collect,1)/randLabelMax
    eigM2CollectMean(n,1:end) = eigM2Collect((n-1)*randLabelMax+1,1:end);
    eigM2CollectMean(n,14) = mean(eigM2Collect((n-1)*randLabelMax+1:n*randLabelMax,14));
    eigM2CollectMean(n,15) = std(eigM2Collect((n-1)*randLabelMax+1:n*randLabelMax,15));
end

figure(111); clf;
hold on;
plot(eigM2CollectMean(:,1),eigM2CollectMean(:,14),'k-')
%set(gca, 'YScale', 'log')
axis([0 1.5 1e-8 5e-2])
hold off;

% NNMlabel = 29;
% NNM(4*(NNMlabel-1)+2,1)



toc
