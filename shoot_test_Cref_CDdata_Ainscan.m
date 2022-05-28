tic

labelAin = 0;

maxZK = 0.51;

for Ain = 0.75:0.05:3
    disp(Ain);
    if Ain <= maxZK
        continue
    else
        
        %     sw = 1;
        %     scd = 1;
        
        wRec = nan(1000,1); % record the omega during the process of calculating the NNMs
        wPer = nan(1000,1);
        
        %     e0 = 1.5; %1.5;
        %     c1 = 0.25;
        %     c2 = scd*0.37;
        %     wChangeFac = 1;
        
        %     d = 0;
        %     d1 = 0.22;
        %     d2 = scd*0.02;
        
        %     f1 = d1;
        %     f2 = d2;
        
        %     if c2 > 0 && d2 > 0
        %         Ac = sqrt(-4*(c2-c1)/3/(d2-d1));
        %     elseif c2 < 0 && d2 < 0
        %         Ac = sqrt(-4*(abs(c2)-c1)/3/(abs(d2)-d1));
        %     end
        %
        %     Ac2 = sqrt(-4*(c2-c1)/3/(d2-d1));
        
        %     etaA = 100;
        %     etaT = 100;
        %     Aed = 10;
        
        %     Ntheo = 10000;
        %     tolTheo = 1e-6;
        %     phiqTheo = zeros(Ntheo,1);
        %     % figure(200); clf;
        %     for nTheo = 1:Ntheo
        %         qTheo = 2*pi/Ntheo*(nTheo-1);
        %         wTheo = sw*sqrt(c1^2+c2^2+2*c1*c2*cos(qTheo));
        %         phiqTheo(nTheo) = asin(c2*sin(qTheo)/wTheo);
        %         if abs(cos(phiqTheo(nTheo))-(c1+c2*cos(qTheo))/wTheo) > tolTheo
        %             phiqTheo(nTheo) = pi-asin(c2*sin(qTheo)/wTheo);
        %         end
        %     end
        % plot(2*pi/Ntheo:2*pi/Ntheo:2*pi,phiqTheo(:)/2/pi)
        % hold off;
        
        
        %     ucNum = 9;
        %     dA = 1/100;
        %     dA2 = dA;
        %     tol1 = 1e-6;
        %     tol2 = 1e-2;
        %     tol3 = 1e-4;
        %     tol4 = 1e-4;
        %     tolA = dA;
        
        
        %     q = 1*2*pi/ucNum;%(floor(ucNum/2))*2*pi/ucNum;
        % sw = +1;
        Delw = sw*sqrt(c1^2+c2^2+2*c1*c2*cos(q));
        phiq0 = asin(c2*sin(q)/Delw);
        if abs(cos(phiq0)-(c1+c2*cos(q))/Delw) > tol1
            phiq0 = pi-asin(c2*sin(q)/Delw);
        end
        w1 = 3*Ain^2/4*(d+d1*cos(phiq0)+d2*cos(phiq0-q));
        phiq1 = -3*Ain^2/8/Delw*(d1*sin(phiq0)+d2*sin(phiq0-q));
        w0 = Delw+e0;
        w = w0+w1;
        Tin = abs(2*pi/w);
        % Tin = abs(2*pi/w0);
        
        %     DOF = 4*ucNum;
        NNM = nan(2,DOF);
        %     timeN = 1000;
        %     critLoopNum = 10000;
        %     nu = 1/4;
        
        Tk = Tin;
        dTk = Tk/timeN;
        
        % initialize the NNM using shooting method
        zk = zeros(timeN+1,DOF);
        gidt = zeros(7,DOF);
        CorrectionList = zeros(DOF,2);
        CorrectionLabel = 0;
        
        maxzk = Ain*ones(DOF,1); % the peak of the disp
        shiftzk = zeros(timeN+1,DOF); % shift the displacement such that the time starts from the maximum position
        shiftLabel = zeros(DOF,1); % record the shifting label
        shiftT = zeros(DOF,1); % record the shifting time
        IntraLabel = zeros(DOF,1); % record the intracell label shift
        IntraPhase = zeros(DOF,1); % record the intracell phase shift
        InterLabel = zeros(DOF,1); % record the intercell label shift
        InterPhase = zeros(DOF,1); % record the intercell phase shift
        
        for n = 1:ucNum
            if sw == -1
                zk(1,4*(n-1)+1) = Ain*sin(+q*n-phiq1); % initial position. here, +pi phase simply means we want to find the lower band, due to the particle hole symmetry!
                zk(1,4*(n-1)+2) = Ain*cos(+q*n-phiq1);
                zk(1,4*(n-1)+3) = Ain*sin(+q*n+phiq0+phiq1);
                zk(1,4*(n-1)+4) = Ain*cos(+q*n+phiq0+phiq1);
            elseif sw == +1
                zk(1,4*(n-1)+1) = Ain*sin(-q*n+phiq0+phiq1); % initial position. here, +pi phase simply means we want to find the lower band, due to the particle hole symmetry!
                zk(1,4*(n-1)+2) = Ain*cos(-q*n+phiq0+phiq1);
                zk(1,4*(n-1)+3) = Ain*sin(-q*n-phiq1);
                zk(1,4*(n-1)+4) = Ain*cos(-q*n-phiq1);
            end
        end
        
        loopA = 0;
        
        while 1
            
            loopA = loopA+1;
            disp(loopA);
            if loopA ~= 1
                if abs((maxzk(1)-Aed)/Aed) < 3*dA
                    break
                else
                    if abs(maxzk(1)) < Ac*0.5
                        zk = zk*(1+5*dA2);
                    else
                        zk = zk*(1+dA2);
                    end
                    if loopA >= 3 % if TkRec has more than 2 records of Tk
                        % wPer = (wRec(loopA-1)-wRec(loopA-2))/wRec(loopA-2); % calculate how much Tk has been changed
                        w = wRec(loopA-1)+(wRec(loopA-1)-wRec(loopA-2))*wChangeFac; % adjust the Tk for the next loop
                        Tk = 2*pi/w;
                        dTk = Tk/timeN;
                    end
                end
            end
            
            loopNum = 0;
            
            while 1
                loopNum = loopNum+1;
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
                
                Hz = zk(timeN+1,:)-zk(1,:);
                errorFunc = abs(sum(abs(Hz)))/DOF;
                
                stdSum = 0;
                for n21 = 1:4
                    stdSum = stdSum+std(InterPhase(n21:4:4*(ucNum-1)+n21));
                    stdSum = stdSum+std(IntraPhase(n21:4:4*(ucNum-1)+n21));
                end
                stdSum = stdSum+std(maxzk);
                
                if stdSum > tol4 && loopNum > critLoopNum
                    break
                else
                    if errorFunc > tol1
                        if mod(loopNum,500) == 0
                            disp(errorFunc)
                        end
                        
                        periodDiff = -transpose(Hz);
                        
                        outM = expm(M2)-eye(DOF,DOF);
                        
                        gf = zeros(DOF,1);
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
                            gf(4*n1-3) = +(+e0*zk(timeN+1,4*n1-2)+d*zk(timeN+1,4*n1-2)^3+c1*zk(timeN+1,4*n1+0)+f1*zk(timeN+1,4*n1+0)^3+c2*zk(timeN+1,4*n1temp-4)+f2*zk(timeN+1,4*n1temp-4)^3);
                            gf(4*n1-2) = -(+e0*zk(timeN+1,4*n1-3)+d*zk(timeN+1,4*n1-3)^3+c1*zk(timeN+1,4*n1-1)+d1*zk(timeN+1,4*n1-1)^3+c2*zk(timeN+1,4*n1temp-5)+d2*zk(timeN+1,4*n1temp-5)^3);
                            gf(4*n1-1) = +(+e0*zk(timeN+1,4*n1-0)+d*zk(timeN+1,4*n1-0)^3+c1*zk(timeN+1,4*n1-2)+f1*zk(timeN+1,4*n1-2)^3+c2*zk(timeN+1,4*n2temp+2)+f2*zk(timeN+1,4*n2temp+2)^3);
                            gf(4*n1-0) = -(+e0*zk(timeN+1,4*n1-1)+d*zk(timeN+1,4*n1-1)^3+c1*zk(timeN+1,4*n1-3)+d1*zk(timeN+1,4*n1-3)^3+c2*zk(timeN+1,4*n2temp+1)+d2*zk(timeN+1,4*n2temp+1)^3);
                        end
                        outM(:,DOF) = gf(:);
                        
                        deltazkT = outM\periodDiff;
                        deltazk = deltazkT(1:DOF-1)/etaA;
                        deltaTk = deltazkT(DOF)/etaT;
                        
                        zk(1,1:DOF-1) = zk(1,1:DOF-1)+transpose(deltazk);
                        Tk = Tk+deltaTk;
                        dTk = Tk/timeN;
                        
                        CorrectionLabel = CorrectionLabel+1;
                        CorrectionList(:,CorrectionLabel) = deltazkT;
                        
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
                                %                     if abs(IntraLabel(4*(n1-1)+n22)-timeN)/timeN < tol2 || abs(IntraLabel(4*(n1-1)+n22))/timeN < tol2
                                %                         IntraLabel(4*(n1-1)+n22) = 0;
                                %                     end
                                IntraPhase(4*(n1-1)+n22) = mod((shiftT(4*(n1-1)+n22)-shiftT(4*(n1-1)+1))/Tk,1);
                                %                     if abs(IntraPhase(4*(n1-1)+n22)-1)/1 < tol3 || abs(IntraLabel(4*(n1-1)+n22))/1 < tol3
                                %                         IntraPhase(4*(n1-1)+n22) = 0;
                                %                     end
                                if n1 == ucNum
                                    n2temp = n1-ucNum;
                                else
                                    n2temp = n1;
                                end
                                InterLabel(4*(n1-1)+n22) = mod((shiftLabel(4*(n2temp-0)+n22)-shiftLabel(4*(n1-1)+n22)),timeN);
                                %                     if abs(InterLabel(4*(n1-1)+n22)-timeN)/timeN < tol2 || abs(InterLabel(4*(n1-1)+n22))/timeN < tol2
                                %                         InterLabel(4*(n1-1)+n22) = 0;
                                %                     end
                                InterPhase(4*(n1-1)+n22) = mod((shiftT(4*(n2temp-0)+n22)-shiftT(4*(n1-1)+n22))/Tk,1);
                                %                     if abs(InterPhase(4*(n1-1)+n22)-1)/1 < tol3 || abs(InterPhase(4*(n1-1)+n22))/1 < tol3
                                %                         InterPhase(4*(n1-1)+n22) = 0;
                                %                     end
                            end
                        end
                        
                    else
                        break
                    end
                end
            end
            if ~isnan(errorFunc)
                
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
                        IntraPhase(4*(n1-1)+n22) = mod((shiftT(4*(n1-1)+n22)-shiftT(4*(n1-1)+1))/Tk,1);
                        if n1 == ucNum
                            n2temp = n1-ucNum;
                        else
                            n2temp = n1;
                        end
                        InterLabel(4*(n1-1)+n22) = mod((shiftLabel(4*(n2temp-0)+n22)-shiftLabel(4*(n1-1)+n22)),timeN);
                        InterPhase(4*(n1-1)+n22) = mod((shiftT(4*(n2temp-0)+n22)-shiftT(4*(n1-1)+n22))/Tk,1);
                    end
                end
                
                NNM(4*(loopA-1)+1,1:DOF) = zk(1,1:DOF);
                NNM(4*(loopA-1)+2,1) = maxzk(1);
                NNM(4*(loopA-1)+2,2) = 2*pi/Tk;
                NNM(4*(loopA-1)+2,3) = errorFunc;
                NNM(4*(loopA-1)+2,4) = sw;
                NNM(4*(loopA-1)+2,5) = scd;
                NNM(4*(loopA-1)+2,6) = ucNum;
                NNM(4*(loopA-1)+2,7) = q;
                NNM(4*(loopA-1)+2,8) = c1;
                NNM(4*(loopA-1)+2,9) = c2;
                NNM(4*(loopA-1)+2,10) = d1;
                NNM(4*(loopA-1)+2,11) = d2;
                NNM(4*(loopA-1)+2,12) = Ac;
                NNM(4*(loopA-1)+3,1:DOF) = transpose(InterPhase(1:DOF));
                NNM(4*(loopA-1)+4,1:DOF) = transpose(IntraPhase(1:DOF));
                wRec(loopA) = 2*pi/Tk;
                if loopA >= 2
                    wPer(loopA-1) = (wRec(loopA)-wRec(loopA-1))/wRec(loopA-1); % calculate how much Tk has been changed
                end
            else
                break
            end
            
            if maxzk(1)^2*max(d1,d2)/e0 > 0.1
                n1 = 1;
                [maxzk1,maxzk1Label] = max(zk(1:timeN+1,4*(n1-1)+1));
                cosCompare = zeros(timeN+1,1);
                for n2 = 1:timeN+1
                    cosCompare(n2) = maxzk1*cos((n2-maxzk1Label)/timeN*2*pi);
                end
                %         figure(loopA); clf;
                %
                %         hold on;
                %         plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+1))
                %         plot(dTk*(1:timeN+1),cosCompare(1:timeN+1))
                %         %     plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+2))
                %         %     plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+3))
                %         %     plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+4))
                %         axis([0 Tk -maxzk(1) maxzk(1)])
                %         hold off;
            end
            
        end
        
        
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
        
        % disp(IntraPhase(3))
        
        if ~isnan(NNM(1,1))
            labelAin = labelAin+1;
            filename1 = strcat('NNMtemp',num2str(labelAin));
            save(filename1,'NNM');
            maxZK = max(maxzk(1),maxZK);
            %         break
        end
    end
end


toc
