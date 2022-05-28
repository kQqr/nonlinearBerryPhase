tic

Aed = 0.6;
loopA = 0;
while 1
    loopA = loopA+1;
    if loopA ~= 1
        if abs((maxzk(1)-Aed)/Aed) < 3*dA
            break
        else
            zk = zk*(1+dA);
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
end


for n1 = 1:ucNum
    figure(n1); clf;
    hold on;
    plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+1))
    plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+2))
    plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+3))
    plot(dTk*(1:timeN+1),zk(1:timeN+1,4*(n1-1)+4))
    axis([0 Tk -maxzk(1) maxzk(1)])
    hold off;
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





toc
