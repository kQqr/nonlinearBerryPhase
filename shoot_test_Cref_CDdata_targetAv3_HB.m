tic

revSign = 1; % if its -1, then flip sign of dA

figSee = 0;
sw = +1;
scd = +1;
for ucNum = 75:2:99
    
    DOF = 4*ucNum;
    tol1 = 1e-6;
    
    timeN = 1000;
    critLoopNum = 3000;
    nu = 1/4;
    wRecSwitch = 0;
    for nq = floor(ucNum/2) % 1:floor(ucNum/2)  % 23:  11:floor(ucNum/2)
        if gcd(ucNum,nq) == 1
            NNM = nan(2,DOF);
            labelNNM2 = 0;
            if sw == 1 && scd == 1
                filename1 = strcat('NNM',num2str(nq),'_',num2str(ucNum),'pp');
            elseif sw == -1 && scd == 1
                filename1 = strcat('NNM',num2str(nq),'_',num2str(ucNum),'mp');
            elseif sw == 1 && scd == -1
                filename1 = strcat('NNM',num2str(nq),'_',num2str(ucNum),'pm');
            elseif sw == -1 && scd == -1
                filename1 = strcat('NNM',num2str(nq),'_',num2str(ucNum),'mm');
            end
            % nq = 2;
            q = nq*2*pi/ucNum;
            % Ain = 0.95;
            for Ain = 1.1 % [0.2 0.4 0.85 0.894427190999916 0.95] % [0.2 0.4 0.85 0.894427190999916 0.95 1.1]
                Atarget = Ain;
                if revSign == -1
                    dA = revSign*1/100;
                else
                    dA = revSign*3/1000;
                end
                
                
                wRec = nan(1000,1); % record the omega during the process of calculating the NNMs
                wPer = nan(1000,1);
                
                e0 = 1.5; %1.5;
                c1 = 0.25;
                c2 = scd*0.37;
                wChangeFac = 1;
                
                d = 0;
                d1 = 0.22;
                d2 = scd*0.02;
                
                f1 = d1;
                f2 = d2;
                
                if c2 > 0 && d2 > 0
                    Ac = sqrt(-4*(c2-c1)/3/(d2-d1));
                elseif c2 < 0 && d2 < 0
                    Ac = sqrt(-4*(abs(c2)-c1)/3/(abs(d2)-d1));
                end
                
                Ac2 = sqrt(-4*(c2-c1)/3/(d2-d1));
                
                
                
                etaA = 300;
                etaT = 10;
                Aed = 10;
                
                
                
                Ntheo = 10000;
                tolTheo = 1e-6;
                phiqTheo = zeros(Ntheo,1);
                for nTheo = 1:Ntheo
                    qTheo = 2*pi/Ntheo*(nTheo-1);
                    wTheo = sw*sqrt(c1^2+c2^2+2*c1*c2*cos(qTheo));
                    phiqTheo(nTheo) = asin(c2*sin(qTheo)/wTheo);
                    if abs(cos(phiqTheo(nTheo))-(c1+c2*cos(qTheo))/wTheo) > tolTheo
                        phiqTheo(nTheo) = pi-asin(c2*sin(qTheo)/wTheo);
                    end
                end
                
                
                tol2 = 1e-2;
                tol3 = 1e-4;
                tol4 = 1e-1;
                
                c1A = c1+3/4*d1*Ain^2;
                c2A = c2+3/4*d2*Ain^2;
                Delw = sw*sqrt(c1A^2+c2A^2+2*c1A*c2A*cos(q));
                phiq0 = asin(c2A*sin(q)/Delw);
                if abs(cos(phiq0)-(c1A+c2A*cos(q))/Delw) > tol1
                    phiq0 = pi-asin(c2A*sin(q)/Delw);
                end
                w1 = 0; %3*Ain^2/4*(d+d1*cos(phiq0)+d2*cos(phiq0-q));
                phiq1 = 0; %-3*Ain^2/8/Delw*(d1*sin(phiq0)+d2*sin(phiq0-q));
                w0 = Delw+e0;
                w = w0+w1;
                Tin = abs(2*pi/w);
                % Tin = abs(2*pi/w0);
                
                
                
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
                        if abs((maxzk(1)-Atarget)/Atarget) < abs(dA)
                            break
                        else
                            if loopA == 2
                                zk = zk*(Atarget/maxzk(1))^(revSign);
                            else
                                if abs((maxzk(1)-Atarget)/Atarget) > 5*abs(dA)
                                    if maxzk(1) > Atarget
                                        zk = zk*(1-3*dA);
                                    elseif maxzk(1) < Atarget
                                        zk = zk*(1+3*dA);
                                    end
                                else
                                    if maxzk(1) > Atarget
                                        zk = zk*(1-dA);
                                    elseif maxzk(1) < Atarget
                                        zk = zk*(1+dA);
                                    end
                                end
                            end
                            
                            if loopA >= 3 % if TkRec has more than 2 records of Tk
                                % wPer = (wRec(loopA-1)-wRec(loopA-2))/wRec(loopA-2); % calculate how much Tk has been changed
                                w = wRec(loopA-1)+wRecSwitch*(wRec(loopA-1)-wRec(loopA-2))*wChangeFac; % adjust the Tk for the next loop
                                Tk = 2*pi/w;
                                dTk = Tk/timeN;
                            end
                        end
                    end
                    
                    loopNum = 0;
                    
                    errorRec = nan(critLoopNum+1,1);
                    
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
                        errorRec(loopNum) = errorFunc;
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
                                
                            else
                                break
                            end
                        end
                    end
                    if ~isnan(errorFunc) && stdSum < tol4
                        
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
                        
                        labelNNM2 = labelNNM2+1;
                        NNM(4*(labelNNM2-1)+1,1:DOF) = zk(1,1:DOF);
                        NNM(4*(labelNNM2-1)+2,1) = maxzk(1);
                        NNM(4*(labelNNM2-1)+2,2) = 2*pi/Tk;
                        NNM(4*(labelNNM2-1)+2,3) = errorFunc;
                        NNM(4*(labelNNM2-1)+2,4) = sw;
                        NNM(4*(labelNNM2-1)+2,5) = scd;
                        NNM(4*(labelNNM2-1)+2,6) = ucNum;
                        NNM(4*(labelNNM2-1)+2,7) = q;
                        NNM(4*(labelNNM2-1)+2,8) = c1;
                        NNM(4*(labelNNM2-1)+2,9) = c2;
                        NNM(4*(labelNNM2-1)+2,10) = d1;
                        NNM(4*(labelNNM2-1)+2,11) = d2;
                        NNM(4*(labelNNM2-1)+2,12) = Ac;
                        NNM(4*(labelNNM2-1)+2,13) = stdSum;
                        NNM(4*(labelNNM2-1)+3,1:DOF) = transpose(InterPhase(1:DOF));
                        NNM(4*(labelNNM2-1)+4,1:DOF) = transpose(IntraPhase(1:DOF));
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
                        
                        if figSee == 1
                            figure(loopA); clf;
                            hold on;
                            plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+1))
                            plot(dTk*(1:timeN+1),cosCompare(1:timeN+1))
                            axis([0 Tk -maxzk(1) maxzk(1)])
                            hold off;
                        end
                    end
                    if stdSum > tol4
                        break
                    end
                    
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
                
            end
            save(filename1,'NNM');
        end
    end
    
    
    % save('NNMtemp','NNM');
    
end

toc
