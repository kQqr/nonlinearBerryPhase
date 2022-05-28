tic
% an example of Gaussian excitation in the topological regime

NtoT = 1;
plotfig10 = 1;
de0 = 0.05;
dc = 0;%0.2;

% if TtoT == 1
    fDrive = 25e-1; % 2.5
    eta = 0.0000;
    % fDrive = 1.5e-2; %1e-2; %5e-4;
    % eta = 0.002;
    fDriveCoeff = 1;
% else
%     fDrive = 7e-2;
%     fDriveCoeff = 1;
% end
if plotfig10 == 1
    pNum = 500; % number of periods 400 for fDrive=2.5
else
    pNum = 500;
end
LineWidth = 3;
blue = [0 0 256]/256;
red = [256 0 0]/256;
yellow = [255 165 0]/256;
gold1 = [255 215 0]/256;
brown = [139 69 19]/256;
alphaValue = 0.4;

ucNum = 120; % 150 (400)
DOF = 4*ucNum;
nu = 1/4;

e0 = 8;
eA = e0*(1+de0);
eB = e0*(1-de0);

w0 = 1;
dw = 0; % if you want dw = e0/20, then need to lower the driving amplitude
d = 0;
c1 = 0.37; %0.25
c2 = 0.25; %0.15
% dc = c2-c1;
d1 = 0.02;
d2 = 0.22;
% dd = d2-d1;

dc1 = dc;
dc2 = dc;
dd1 = dc;
dd2 = dc;

c1n = c1*ones(ucNum,1)+c1*dc1*(rand(ucNum,1)-0.5);
c2n = c2*ones(ucNum,1)+c2*dc2*(rand(ucNum,1)-0.5);
d1n = d1*ones(ucNum,1)+d1*dd1*(rand(ucNum,1)-0.5);
d2n = d2*ones(ucNum,1)+d2*dd2*(rand(ucNum,1)-0.5);

% save('OBCdisorderc10_2','c1n');
% save('OBCdisorderc20_2','c2n');
% save('OBCdisorderd10_2','d1n');
% save('OBCdisorderd20_2','d2n');

Ac = sqrt(-4*(c2-c1)/3/(d2-d1));

% et = e0+abs(dc);
% eb = e0-abs(dc);



f = zeros(DOF,1);
% in the topo regime, no matter we excite x1, y1, or we excite XN, YN, we always get topo excitation! in the non-topo regime, no matter we excite
% x1, y1, or we excite XN, YN, we always get a small excitation, although also exponential increase as (c2/c1)^{-1} (when c1>c2 non-topo)



timeN = pNum*100; % 200000


w = e0+dw;
dt = 2*pi/w/(timeN/pNum);
nT = floor(2*pi/w/dt);

tnghw = 10*nT;% Gaussian half width of the driving signal, in units of the time count, nt
tnStart = 25*nT; % starting time of the shaking force of the Gaussian signal

% timeN = floor(600*nT);


tInterval = timeN/2:timeN; % 10*nT:timeN;
%tInterval = timeN/4+10*nT:timeN/4+timeN/2+10*nT-1;

% eta = 0.000;
% fDrive = 0.085; % bulk + edge modes 
% tInterval = timeN-60*nT:timeN;



inputSite1 = 1;
inputSite2 = 2;
% inputSite1 = DOF-1;
% inputSite2 = DOF-0;
f(inputSite1) = fDrive*fDriveCoeff;
f(inputSite2) = fDrive*fDriveCoeff;
phase = zeros(DOF,1);
phase(inputSite1) = 0;
phase(inputSite2) = pi/2;

z = zeros(timeN,DOF);
gidt = zeros(7,DOF);

fRec = nan(timeN,2); % record the time-domain driving force

for n = 1:timeN
    yi = zeros(7,DOF);
    gFactor = exp(-((n-tnStart)*dt)^2/(tnghw*dt)^2); % the Gaussian shape function
    fRec(n,1) = f(1)*sin(w*n*dt+phase(1))*gFactor;
    fRec(n,2) = f(2)*sin(w*n*dt+phase(2))*gFactor;
    for n20 = 1:7
        if n20 == 1
            yi(1,:) = z(n,:);
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
            
            Fa  = eA*yi(n20,4*n1-2)+ d*yi(n20,4*n1-2)^3;
            F1a = c1n(n1)*yi(n20,4*n1+0)+d1n(n1)*yi(n20,4*n1+0)^3;
            if n1 == 1
                F2a = 0;
            else
                F2a = c2n(n1-1)*yi(n20,4*n1-4)+d2n(n1-1)*yi(n20,4*n1-4)^3;
            end
            Ga  = +Fa+F1a+F2a-f(4*n1-3)*sin(w*n*dt+phase(4*n1-3))*gFactor-eta*w0*yi(n20,4*n1-3);
            
            Fb  = eA*yi(n20,4*n1-3)+ d*yi(n20,4*n1-3)^3;
            F1b = c1n(n1)*yi(n20,4*n1-1)+d1n(n1)*yi(n20,4*n1-1)^3;
            if n1 == 1
                F2b = 0;
            else
                F2b = c2n(n1-1)*yi(n20,4*n1-5)+d2n(n1-1)*yi(n20,4*n1-5)^3;
            end
            Gb  = -Fb-F1b-F2b-f(4*n1-2)*sin(w*n*dt+phase(4*n1-2))*gFactor-eta*w0*yi(n20,4*n1-2);
            
            Fc  = eB*yi(n20,4*n1+0)+ d*yi(n20,4*n1+0)^3;
            F1c = c1n(n1)*yi(n20,4*n1-2)+d1n(n1)*yi(n20,4*n1-2)^3;
            if n1 == ucNum
                F2c = 0;
            else
                F2c = c2n(n1)*yi(n20,4*n1+2)+d2n(n1)*yi(n20,4*n1+2)^3;
            end
            Gc  = +Fc+F1c+F2c-f(4*n1-1)*sin(w*n*dt+phase(4*n1-1))*gFactor-eta*w0*yi(n20,4*n1-1);
            
            Fd  = eB*yi(n20,4*n1-1)+ d*yi(n20,4*n1-1)^3;
            F1d = c1n(n1)*yi(n20,4*n1-3)+d1n(n1)*yi(n20,4*n1-3)^3;
            if n1 == ucNum
                F2d = 0;
            else
                F2d = c2n(n1)*yi(n20,4*n1+1)+d2n(n1)*yi(n20,4*n1+1)^3;
            end
            Gd  = -Fd-F1d-F2d-f(4*n1-0)*sin(w*n*dt+phase(4*n1-0))*gFactor-eta*w0*yi(n20,4*n1-0);

%             gidt(n20,4*n1-3) = +dt*(Ga-eta*Gb);
%             gidt(n20,4*n1-2) = -dt*(Gb+eta*Ga);
%             gidt(n20,4*n1-1) = +dt*(Gc-eta*Gd);
%             gidt(n20,4*n1-0) = -dt*(Gd+eta*Gc);

            gidt(n20,4*n1-3) = dt*Ga;
            gidt(n20,4*n1-2) = dt*Gb;
            gidt(n20,4*n1-1) = dt*Gc;
            gidt(n20,4*n1-0) = dt*Gd;
        end
    end
    
    z(n+1,:) = z(n,:)+(9*gidt(1,:)+64*gidt(3,:)+49*gidt(5,:)+49*gidt(6,:)+9*gidt(7,:))/180;
    
end

figure(6); clf;
hold on;
H1 = plot((1:size(z,1))*dt/(2*pi/e0),z(:,1),'color',red,'linewidth',LineWidth/3);
H2 = plot((1:size(z,1))*dt/(2*pi/e0),z(:,5),'color',blue,'linewidth',LineWidth/3);
H3 = plot((1:size(z,1))*dt/(2*pi/e0),z(:,9),'color',yellow,'linewidth',LineWidth/3);
hold off;
legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('$\Psi_n^{(1)}(t)$','Interpreter','latex')
xticks([0 200 400])
yticks([-20 0 20])
%yticklabels({'-A_c','-0.5','0','0.5','A_c'})
axis([0 400 -25 25])







figure(107); clf;
hold on;
H1 = plot((1:size(z,1))*dt/(2*pi/e0),z(1:size(z,1),4*(10-1)+1),'color',red,'linewidth',LineWidth/3);
H2 = plot((1:size(z,1))*dt/(2*pi/e0),z(1:size(z,1),4*(20-1)+1),'color',blue,'linewidth',LineWidth/3);
H3 = plot((1:size(z,1))*dt/(2*pi/e0),z(1:size(z,1),4*(37-1)+1),'color',yellow,'linewidth',LineWidth/3);
hold off;
legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('$\Psi_n^{(1)}(t)$','Interpreter','latex')







figure(7); clf;
hold on;
H1 = plot((1:size(z,1))*dt/(2*pi/e0),z(1:size(z,1),4*(1-1)+3),'color',red,'linewidth',LineWidth/3);
% H2 = plot((1:size(z,1))*dt/(2*pi/e0),z(1:size(z,1),4*(20-1)+3),'color',blue,'linewidth',LineWidth/3);
% H3 = plot((1:size(z,1))*dt/(2*pi/e0),z(1:size(z,1),4*(37-1)+3),'color',yellow,'linewidth',LineWidth/3);
hold off;
legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('$\Psi_n^{(1)}(t)$','Interpreter','latex')










figure(14); clf;
hold on;

% [maxz1,maxz1Label] = max(z(:,1));
% snapShotTime = maxz1Label/(2*pi/e0);

% ampz = zeros(timeN+1,2*ucNum);
% for n1 = 1:timeN+1
%     for n = 1:2*ucNum
%         ampz(n1,n) = sqrt((z(n1,2*(n-1)+1))^2+(z(n1,2*(n-1)+2))^2);
%     end
% end

ampzMax = zeros(2*ucNum,1);
for n = 1:2*ucNum
    % ampzMax(n) = max(ampz(:,n));
    ampzMax(n) = max(z(:,2*(n-1)+1));
end
    
% plot(1:ucNum,z(maxz1Label,1:4:4*(ucNum-1)+1))
% plot(1:ucNum,z(maxz1Label,2:4:4*(ucNum-1)+2))
% plot(1:ucNum,z(maxz1Label,3:4:4*(ucNum-1)+3))
% plot(1:ucNum,z(maxz1Label,4:4:4*(ucNum-1)+4))

for n = 1:ucNum
    H1 = plot(n,ampzMax(2*n-1),'ro','markersize',5,'linewidth',2);
    %plot([n n],[0 ampzMax(2*n-1)],'r-','linewidth',2)
    H2 = plot(n,ampzMax(2*n-0),'bo','markersize',5,'linewidth',2);
    %plot([n n],[0 ampzMax(2*n-0)],'b-','linewidth',2)
end
legend([H1 H2],{'${\rm max\,}({\rm Re\,}\Psi_n^{(1)})$','${\rm max\,}({\rm Re\,}\Psi_n^{(2)})$'},'Location','northeast','FontSize',16,'Interpreter','latex');
% axis([0 ucNum 0 1.1*max(ampzMax)])
axis([0 ucNum 1e-1 1.1*max(ampzMax)])
xlabel('Site Labeling','Interpreter','latex')
% ylabel('${\rm max\,}|\Psi_n^{(1)}|$','Interpreter','latex')
% xticks([0 15 30 45])
% yticks([0 Ac/2 Ac])
% yticklabels({'0','A_c/2','A_c'})
hold off;
set(gca,'linewidth',2)
set(gca,'FontSize',19)
set(gca, 'YScale', 'log')



% ampz = zeros(DOF,1);
% for n1 = 1:ucNum
%     ampz(4*n1-3) = (max(z(tInterval,4*n1-3))-min(z(tInterval,4*n1-3)))/2;
%     ampz(4*n1-2) = (max(z(tInterval,4*n1-2))-min(z(tInterval,4*n1-2)))/2;
%     ampz(4*n1-1) = (max(z(tInterval,4*n1-1))-min(z(tInterval,4*n1-1)))/2;
%     ampz(4*n1-0) = (max(z(tInterval,4*n1-0))-min(z(tInterval,4*n1-0)))/2;
% end
% figure(7); clf;
% hold on;
% plot(1:ucNum,ampz(1:4:DOF-3))
% plot(1:ucNum,ampz(2:4:DOF-2))
% plot(1:ucNum,ampz(3:4:DOF-1))
% plot(1:ucNum,ampz(4:4:DOF-0))
% hold off;
% set(gca, 'YScale', 'log')
% DtestPosi1 = 1;
% DtestPosi2 = floor(ucNum/3);
% DtestPosiDiff = 5;
% DrateDiff1 = (log(ampz(4*DtestPosi1-3)/ampz(4*(DtestPosi1+DtestPosiDiff)-3))/DtestPosiDiff - log(c2/c1))/log(c2/c1);
% disp(DrateDiff1);
% DrateDiff2 = (log(ampz(4*DtestPosi2-3)/ampz(4*(DtestPosi2+DtestPosiDiff)-3))/DtestPosiDiff - log(c2/c1))/log(c2/c1);
% disp(DrateDiff2);
% disp(max(ampz));







figure(8); clf;
hold on;
H1 = plot((1:timeN)*dt/(2*pi/e0),fRec(:,1),'b-','linewidth',LineWidth/2);
% H2 = plot((1:timeN)*dt/(2*pi/e0),fRec(:,2),'r-','linewidth',LineWidth/2);
hold off;
legend([H1],{'${\rm Re\,}S_1(t)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('${\rm Re\,}S_1(t)$','Interpreter','latex')
xticks([0 25 50])
yticks([-fDrive 0 fDrive])
axis([0 2*tnStart/nT -fDrive*1.1 fDrive*1.1])





FreqMeasure = 2*pi/dt/size(tInterval,2)*(0:(size(tInterval,2)/2)); % the ruler to measure the spectrum: in unit of freqUnit
for n1 = 1:size(FreqMeasure,2)
    if FreqMeasure(n1) > e0
        break
    end
end
freqMeasureLabel = 4*n1;

WF = nan(timeN+1,ucNum); % the wave function in the complex form
for n = 1:timeN
    for n1 = 1:ucNum 
        WF(n,n1) = z(n,2*(n1-1)+1)+1i*z(n,2*(n1-1)+2);
    end
end

ampwz = nan(freqMeasureLabel,4*ucNum);% record the FFT but limit to 4*f(e0) cuz we want the 3w mode!

for n1 = 1:4*ucNum
    
    tol6 = 1e-3;
    P1x = 0;
    Q1x = 0;
    P2x = 0;
    Q2x = 0;
    
    P2x = 2/size(tInterval,2)*abs(fft(z(tInterval,n1)));
    P1x = P2x;
    P1x = P2x(1:floor(size(tInterval,2)/2)+1);
    P1x(2:end-1) = P1x(2:end-1);
    
    Q2x = 2/size(tInterval,2)*abs(fft(z(tInterval,n1)));
    Q1x = Q2x;
    Q1x = Q2x(1:floor(size(tInterval,2)/2)+1);
    Q1x(2:end-1) = Q1x(2:end-1);
    
    [maxP1x,maxP1xLabel] = max(P1x);
    
    ampwz(:,n1) = P1x(1:freqMeasureLabel);

end




ampWF = nan(freqMeasureLabel,ucNum);% record the FFT but limit to 4*f(e0) cuz we want the 3w mode!

for n1 = 1:ucNum
    
    tol6 = 1e-3;
    P1x = 0;
    Q1x = 0;
    P2x = 0;
    Q2x = 0;
    
    P2x = 2/size(tInterval,2)*abs(fft(WF(tInterval,n1)));
    P1x = P2x;
    P1x = P2x(1:floor(size(tInterval,2)/2)+1);
    P1x(2:end-1) = P1x(2:end-1);
    
    Q2x = 2/size(tInterval,2)*abs(fft(WF(tInterval,n1)));
    Q1x = Q2x;
    Q1x = Q2x(1:floor(size(tInterval,2)/2)+1);
    Q1x(2:end-1) = Q1x(2:end-1);
    
    [maxP1x,maxP1xLabel] = max(P1x);
    
    ampWF(:,n1) = P1x(1:freqMeasureLabel);

end




% ampwzAbs = nan(size(ampwz,1),2*ucNum); % sqrt(re^2+im^2)
% for n1 = 1:size(ampwz,1)
%     for n2 = 1:2*ucNum
%         ampwzAbs(n1,n2) = sqrt(ampwz(n1,2*(n2-1)+1)^2+ampwz(n1,2*(n2-1)+2)^2);
%     end
% end

Aedge = max(ampwz(1:end,1));


[Bedge2,maxBLabel] = max(ampwz(1:end,4*(1-1)+3));

wtol = 1e-6;
c0 = e0-(c2-c1);
d0 = d-(d2-d1);
x0 = -c0/d0;
wNonlinear = nan(1,3);
label = 0;
for aUp = Aedge
    R2 = (aUp^2-x0)^2+x0^2;
    
    label = label+1;
    fun = @(a) 1./(sqrt(R2-(a.^2-x0).^2).*sqrt(x0+sqrt(R2-(a.^2-x0).^2)));
    intF = integral(fun,0,aUp);
    
    wNonlinear(label,1) = abs(d0)*pi/2/abs(intF);
    wNonlinear(label,2) = sqrt(R2);
    wNonlinear(label,3) = aUp;
    
    if abs(wNonlinear(label,1)) <= wtol
        wNonlinear(label,1) = nan;
    end
end

C0 = e0+(c2-c1);
D0 = d+(d2-d1);
X0 = -C0/D0;
WNonlinear = nan(1,2);
label = 0;
for aUp = Aedge
    R2 = (aUp^2-X0)^2+X0^2;
    label = label+1;
    fun = @(a) 1./(sqrt(R2-(a.^2-X0).^2).*sqrt(X0-sqrt(R2-(a.^2-X0).^2)));
    intF = integral(fun,0,aUp);
    
    WNonlinear(label,1) = abs(D0)*pi/2/abs(intF);
    WNonlinear(label,2) = sqrt(R2);
    WNonlinear(label,3) = aUp;
    
    if abs(WNonlinear(label,1)) <= wtol
        WNonlinear(label,1) = nan;
    end
end



















tol6 = 1e-3;
P1x = 0;
Q1x = 0;
P2x = 0;
Q2x = 0;

FreqMeasure2 = 2*pi/dt/timeN*(0:(timeN/2)); % the ruler to measure the spectrum: in unit of freqUnit

P2x = 1/(tnghw*sqrt(pi)/dt)*abs(fft(fRec(:,1)));
P1x = P2x;
P1x = P2x(1:floor(timeN/2)+1);
P1x(2:end-1) = P1x(2:end-1);

Q2x = 1/(tnghw*sqrt(pi)/dt)*abs(fft(fRec(:,2)));
Q1x = Q2x;
Q1x = Q2x(1:floor(timeN/2)+1);
Q1x(2:end-1) = Q1x(2:end-1);

[maxP1x,maxP1xLabel] = max(P1x);
fftLabelLimit = floor(maxP1xLabel*2);

figure(9); clf;
hold on;

bandgap = sqrt((eA-eB)^2/4+(c1-c2)^2);
bsBot = e0-bandgap;
bsTop = e0+bandgap;
% X = [wNonlinear(1,1); WNonlinear(1,1); WNonlinear(1,1); wNonlinear(1,1)];
X = [bsBot; bsTop; bsTop; bsBot];
Y = [0; 0; 1; 1];
h = fill(X,Y,gold1);
set(h,'EdgeColor','none')
alpha(alphaValue)

H1 = plot(FreqMeasure2(1:fftLabelLimit),P1x(1:fftLabelLimit)/max(P1x(1:fftLabelLimit)),'b-','linewidth',LineWidth,'color',brown);
H2 = plot(FreqMeasure(1:fftLabelLimit),ampwz(1:fftLabelLimit,1)/max(ampwz(1:fftLabelLimit,1)),'b-','linewidth',3);
% plot([+c1-c2 +c1-c2],[0 3e-4],'k--','linewidth',LineWidth)
% plot([-c1+c2 -c1+c2],[0 3e-4],'k--','linewidth',LineWidth)
hold off;
% legend([H1 H2],{'$S_1(\omega)/\max(S_1)$','$\Psi_1^{(1)}(\omega)/\max(\Psi_1^{(1)})$'},'Location','northeast','FontSize',16,'Interpreter','latex');
legend([H1 H2],{'${\rm Re\,}S_1(\omega)$','${\rm Re\,}\Psi_1^{(1)}(\omega)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$\omega$','Interpreter','latex')
ylabel('FFT','Interpreter','latex')
% xticks([e0-abs(c1-c2)*10 e0 e0+abs(c1-c2)*10])
xticks([e0-1 e0 e0+1])
yticks([0 0.5 1])
axis([e0-1 e0+1 0 1])
% set(gca, 'YScale', 'log')






tol6 = 1e-3;
P1x = 0;
Q1x = 0;
P2x = 0;
Q2x = 0;
FreqMeasure2 = 2*pi/dt/timeN*(0:(timeN/2)); % the ruler to measure the spectrum: in unit of freqUnit
P2x = 1/(tnghw*sqrt(pi)/dt)*abs(fft(fRec(:,1)));
P1x = P2x;
P1x = P2x(1:floor(timeN/2)+1);
P1x(2:end-1) = P1x(2:end-1);
Q2x = 1/(tnghw*sqrt(pi)/dt)*abs(fft(fRec(:,2)));
Q1x = Q2x;
Q1x = Q2x(1:floor(timeN/2)+1);
Q1x(2:end-1) = Q1x(2:end-1);
[maxP1x,maxP1xLabel] = max(P1x);
fftLabelLimit = floor(maxP1xLabel*2);
fig91yLim = FreqMeasure(fftLabelLimit);
ucNumTot = 45;
ucNumsc = ucNumTot*2;
imagescFFTz = zeros(fftLabelLimit,4);
for ntemp3 = 1:ucNumsc
    imagescFFTz(:,ntemp3) = ampwz(1:fftLabelLimit,2*ntemp3-1);
end

figure(91); clf;
hold on;

bandgap = sqrt((eA-eB)^2/4+(c1-c2)^2);
bsBot = e0-bandgap;
bsTop = e0+bandgap;
% X = [wNonlinear(1,1); WNonlinear(1,1); WNonlinear(1,1); wNonlinear(1,1)];
%imagesc(1-0.5:1:ucNumTot-0.5,FreqMeasure(1:fftLabelLimit),imagescFFTz(:,1:2:end),[0 5])
imagesc(1-0.5:1:ucNumTot-0.5,FreqMeasure(1:fftLabelLimit),imagescFFTz(:,1:2:end),[0 max(max(imagescFFTz(:,1:2:end)))])
plot([0 ucNumTot],[bsBot bsBot],'w--','linewidth',1)
plot([0 ucNumTot],[bsTop bsTop],'w--','linewidth',1)

axis([0 ucNumTot 7.5 8.5])
xlabel('$n$','Interpreter','latex')
% yticks([-Ac -0.5 0 0.5 Ac])
% yticklabels({'-A_c','-0.5','0','0.5','A_c'})
xticks([15-0.5 30-0.5 45-0.5])
xticklabels([15 30 45])
ylabel('$\omega$','Interpreter','latex')
yticks([7.5 8 8.5])
hold off;
colorbar
colormap(hot)
%colorbar('Ticks',[-5,-2,1,4,7],...
%'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
colorbar('Ticks',[0 5])
% legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
% set(gca,'linewidth',2)
% set(gca,'FontSize',19)





figure(109); clf;
hold on;

bandgap = sqrt((eA-eB)^2/4+(c1-c2)^2);
bsBot = e0-bandgap;
bsTop = e0+bandgap;
% X = [wNonlinear(1,1); WNonlinear(1,1); WNonlinear(1,1); wNonlinear(1,1)];
X = [bsBot; bsTop; bsTop; bsBot];
Y = [0; 0; 1; 1];
h = fill(X,Y,gold1);
set(h,'EdgeColor','none')
alpha(alphaValue)

%H1 = plot(FreqMeasure2(1:fftLabelLimit),P1x(1:fftLabelLimit)/max(P1x(1:fftLabelLimit)),'b-','linewidth',LineWidth,'color',brown);
H2 = plot(FreqMeasure(1:fftLabelLimit),ampwz(1:fftLabelLimit,30)/max(ampwz(1:fftLabelLimit,30)),'b-','linewidth',3);
% plot([+c1-c2 +c1-c2],[0 3e-4],'k--','linewidth',LineWidth)
% plot([-c1+c2 -c1+c2],[0 3e-4],'k--','linewidth',LineWidth)
hold off;
% legend([H1 H2],{'$S_1(\omega)/\max(S_1)$','$\Psi_1^{(1)}(\omega)/\max(\Psi_1^{(1)})$'},'Location','northeast','FontSize',16,'Interpreter','latex');
%legend([H1 H2],{'${\rm Re\,}S_1(\omega)$','${\rm Re\,}\Psi_1^{(1)}(\omega)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$\omega$','Interpreter','latex')
ylabel('FFT','Interpreter','latex')
xticks([e0-abs(c1-c2)*2 e0 e0+abs(c1-c2)*2])
yticks([0 0.5 1])
axis([e0-abs(c1-c2)*2 e0+abs(c1-c2)*2 0 1])
% set(gca, 'YScale', 'log')





[e0DiffMin,e0DiffMinLabel] = min(abs(FreqMeasure-e0*ones(size(FreqMeasure,1),size(FreqMeasure,2))));
maxALabel = e0DiffMinLabel;
% [maxA,maxALabel] = max(ampwz(:,1));
Bedge = ampwz(maxALabel,3);
% FreqMeasure(maxALabel)
figure(10); clf;
hold on;

% below we look for theoretical fitting
tol = 1e-6;
zeta = (-1+1i*sqrt(3))/2;
topoAmpA = zeros(ucNum,1);
topoAmpA(1) = Aedge;
topoAmpB = zeros(ucNum,1);
topoAmpB(1) = Bedge;
ACube = 1;
BCube = 0;
CCube = 4*c2/3/d2;
CCubeB = 4*c1/3/d1;

for n = 1:ucNum-1
    Apre = topoAmpA(n);
    DCube = -d1/d2*Apre^3-4*c1/3/d2*Apre;
    Delta0 = BCube^2-3*ACube*CCube;
    Delta1 = 2*BCube^3-9*ACube*BCube*CCube+27*ACube^2*DCube;
    Ctemp = ((Delta1+sqrt(Delta1^2-4*Delta0^3))/2)^(1/3);
    x = zeros(3,1);
    for n1 = 1:3
        x(n1) = -1/3/ACube*(BCube+zeta^n1*Ctemp+Delta0/Ctemp/zeta^n1);
        if abs(imag(x(n1))) < tol
            topoAmpA(n+1) = real(x(n1));
        end
    end
    %     if abs(d2*Apre^2/c2) < tol
    %         topoAmpA(n+1) = Apre*c1/c2;
    %     end
    
    % topoAmpB(n+1) = topoAmpB(n)*(c2+3*d2/4*Apre^2)/(c1+3*d1/4*Apre^2);
    %     Bpre = topoAmpB(n);
    %     DCube = -d2/d1*Bpre^3-4*c2/3/d1*Bpre;
    %     Delta0 = BCube^2-3*ACube*CCubeB;
    %     Delta1 = 2*BCube^3-9*ACube*BCube*CCubeB+27*ACube^2*DCube;
    %     Ctemp = ((Delta1+sqrt(Delta1^2-4*Delta0^3))/2)^(1/3);
    %     x = zeros(3,1);
    %     for n1 = 1:3
    %         x(n1) = -1/3/ACube*(BCube+zeta^n1*Ctemp+Delta0/Ctemp/zeta^n1);
    %         if abs(imag(x(n1))) < tol && real(x(n1)) > Bpre*c2/c1
    %             topoAmpB(n+1) = real(x(n1));
    %         end
    %     end
    
%     Bpre = topoAmpB(n);
%     DCube = -d1/d2*Bpre^3-4*c1/3/d2*Bpre;
%     Delta0 = BCube^2-3*ACube*CCube;
%     Delta1 = 2*BCube^3-9*ACube*BCube*CCube+27*ACube^2*DCube;
%     Ctemp = ((Delta1+sqrt(Delta1^2-4*Delta0^3))/2)^(1/3);
%     x = zeros(3,1);
%     for n1 = 1:3
%         x(n1) = -1/3/ACube*(BCube+zeta^n1*Ctemp+Delta0/Ctemp/zeta^n1);
%         if abs(imag(x(n1))) < tol
%             topoAmpB(n+1) = real(x(n1));
%         end
%     end
%     %     if abs(d2*Apre^2/c2) < tol
%     %         topoAmpA(n+1) = Apre*c1/c2;
%     %     end
    
end


H1 = plot(1:ucNum,ampwz(maxALabel,4*(1-1)+1:4:4*(ucNum-1)+1),'r-','linewidth',LineWidth);
% H1 = plot(1:ucNum,ampwzAbs(maxALabel,2*(1-1)+1:2:2*(ucNum-1)+1),'r-','linewidth',LineWidth);
% plot(1:ucNum,ampwz(maxALabel,4*(1-1)+2:4:4*(ucNum-1)+2),'b-')
H2 = plot(1:ucNum,ampwz(maxALabel,4*(1-1)+3:4:4*(ucNum-1)+3),'b-','linewidth',LineWidth);
% H2 = plot(1:ucNum,ampwzAbs(maxALabel,2*(1-1)+2:2:2*(ucNum-1)+2),'b-','linewidth',LineWidth);

if de0 == 0
    H3 = plot(1:ucNum,topoAmpA(:),'k--','linewidth',LineWidth/2);
end
% H4 = plot(1:ucNum,topoAmpB(:),'k--','linewidth',LineWidth/2);
hold off;
if de0 == 0
    legend([H1 H3 H2],{'${\rm Re\,}\Psi^{(1)}_n(\epsilon_0)$','$\psi^{(1)}_n(\epsilon_0)$','${\rm Re\,}\Psi_n^{(2)}(\epsilon_0)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
else
    legend([H1 H2],{'${\rm Re\,}\Psi^{(1)}_n(\epsilon_0)$','${\rm Re\,}\Psi_n^{(2)}(\epsilon_0)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
end
% legend([H1 H2],{'$|\Psi^{(1)}_n(\epsilon_0)|$','$|\Psi_n^{(2)}(\epsilon_0)|$'},'Location','northeast','FontSize',16,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('Site labeling','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')
% xticks([0 ucNum/3 2*ucNum/3 ucNum])

% yticklabels({'10^{-6}','10^{-3}','A_c'})
% axis([1 ucNum 1e-2 20])
axis([1 45 5e-5 20])
xticks([15 30 45])
yticks([1e-5 1e-3 1e-1 10])
set(gca, 'YScale', 'log')






% tInterval0 = 5*nT:timeN;
% FreqMeasure = 2*pi/dt/size(tInterval0,2)*(0:(size(tInterval0,2)/2)); % the ruler to measure the spectrum: in unit of freqUnit
% tol6 = 1e-3;
% P1x = 0;
% Q1x = 0;
% P2x = 0;
% Q2x = 0;
% 
% P2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,10)));
% P1x = P2x;
% P1x = P2x(1:floor(size(tInterval0,2)/2)+1);
% P1x(2:end-1) = P1x(2:end-1);
% 
% Q2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,10)));
% Q1x = Q2x;
% Q1x = Q2x(1:floor(size(tInterval0,2)/2)+1);
% Q1x(2:end-1) = Q1x(2:end-1);
% 
% [maxP1x,maxP1xLabel] = max(P1x);
% 
% figure(11); clf;
% hold on;
% plot(FreqMeasure(1:maxP1xLabel*10),P1x(1:maxP1xLabel*10),'k-','linewidth',3)
% hold off;








% fftSite = 4*(91-1)+1;
fftSite = 4*(1-1)+9;
tInterval0 = timeN/2:timeN; % 5*nT:timeN;
FreqMeasure = 2*pi/dt/size(tInterval0,2)*(0:(size(tInterval0,2)/2)); % the ruler to measure the spectrum: in unit of freqUnit
tol6 = 1e-3;
P1x = 0;
Q1x = 0;
P2x = 0;
Q2x = 0;

P2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,fftSite)));
P1x = P2x;
P1x = P2x(1:floor(size(tInterval0,2)/2)+1);
P1x(2:end-1) = P1x(2:end-1);

Q2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,fftSite)));
Q1x = Q2x;
Q1x = Q2x(1:floor(size(tInterval0,2)/2)+1);
Q1x(2:end-1) = Q1x(2:end-1);

[maxP1x,maxP1xLabel] = max(P1x);

figure(12); clf;
hold on;
plot(FreqMeasure(1:min(maxP1xLabel*10,size(FreqMeasure,2))),P1x(1:min(maxP1xLabel*10,size(FreqMeasure,2))),'k-','linewidth',3)
hold off;










tInterval0 = timeN/2:timeN;%5*nT:timeN;
FreqMeasure = 2*pi/dt/size(tInterval0,2)*(0:(size(tInterval0,2)/2)); % the ruler to measure the spectrum: in unit of freqUnit
tol6 = 1e-3;
P1x = 0;
Q1x = 0;
P2x = 0;
Q2x = 0;

P2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,3)));
P1x = P2x;
P1x = P2x(1:floor(size(tInterval0,2)/2)+1);
P1x(2:end-1) = P1x(2:end-1);

Q2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,3)));
Q1x = Q2x;
Q1x = Q2x(1:floor(size(tInterval0,2)/2)+1);
Q1x(2:end-1) = Q1x(2:end-1);

[maxP1x,maxP1xLabel] = max(P1x);

figure(13); clf;
hold on;
plot(FreqMeasure(1:min(maxP1xLabel*10,size(FreqMeasure,2))),P1x(1:min(maxP1xLabel*10,size(FreqMeasure,2))),'k-','linewidth',3)
hold off;







% (c*atemp1+g*atemp1^3)/(eta*e0*atemp1)
% (c*atemp1+g*atemp1^3)/(eta*e0*Atemp1)

figure(10001); clf;
hold on;

t1Label = floor(timeN/10*9);
t2Label = t1Label+450;
tRange = t1Label:t2Label;
ucNum1 = ucNum/2;

H1 = plot(tRange*dt/(2*pi/e0),z(tRange,4*(ucNum1-1)+1),'color',red,'linewidth',LineWidth/3);
[max1,maxLabel1] = max(z(tRange,4*(ucNum1-1)+1));
phase1 = mod(maxLabel1*dt/(2*pi/e0),1);

H2 = plot(tRange*dt/(2*pi/e0),z(tRange,4*(ucNum1-1)+2),'color',red,'linewidth',LineWidth/3);
[max2,maxLabel2] = max(z(tRange,4*(ucNum1-1)+2));
phase2 = mod(maxLabel2*dt/(2*pi/e0),1);

H3 = plot(tRange*dt/(2*pi/e0),z(tRange,4*(ucNum1-1)+3),'color',blue,'linewidth',LineWidth/3);
[max3,maxLabel3] = max(z(tRange,4*(ucNum1-1)+3));
phase3 = mod(maxLabel3*dt/(2*pi/e0),1);

H4 = plot(tRange*dt/(2*pi/e0),z(tRange,4*(ucNum1-1)+4),'color',yellow,'linewidth',LineWidth/3);
[max4,maxLabel4] = max(z(tRange,4*(ucNum1-1)+4));
phase4 = mod(maxLabel4*dt/(2*pi/e0),1);

H5 = plot(tRange*dt/(2*pi/e0),z(tRange,4*(ucNum1-1)+5),'color',red,'linewidth',LineWidth/3);
[max1,maxLabel1] = max(z(tRange,4*(ucNum1-1)+1));
phase1 = mod(maxLabel1*dt/(2*pi/e0),1);

[max5,maxLabel5] = max(z(tRange,4*(ucNum1-1)+5));
phase5 = mod(maxLabel5*dt/(2*pi/e0),1);

[max7,maxLabel7] = max(z(tRange,4*(ucNum1-1)+7));
phase7 = mod(maxLabel7*dt/(2*pi/e0),1);

[max9,maxLabel9] = max(z(tRange,4*(ucNum1-1)+9));
phase9 = mod(maxLabel9*dt/(2*pi/e0),1);

[max11,maxLabel11] = max(z(tRange,4*(ucNum1-1)+11));
phase11 = mod(maxLabel11*dt/(2*pi/e0),1);

[max15,maxLabel15] = max(z(tRange,4*(ucNum1-1)+15));
phase15 = mod(maxLabel15*dt/(2*pi/e0),1);

hold off;
% legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('$\Psi_n^{(1)}(t)$','Interpreter','latex')

Phase2 = phase2-phase1;
Phase3 = phase3-phase1;
Phase4 = phase4-phase1;
Phase5 = phase5-phase1;
Phase7 = phase7-phase1;
Phase9 = phase9-phase1;
Phase11 = phase11-phase1;
Phase15 = phase15-phase1;

c = c2-c1;
d = d2-d1;
g = d*3/4;

% (c*atemp1+g*atemp1^3-eta*e0*Atemp1)/(eta*e0*Atemp1)
% % (c*atemp2+g*atemp2^3-eta*e0*Atemp1)
% (c*Atemp1+g*Atemp1^3+eta*e0*atemp1)/(eta*e0*atemp1)
% % (c*Atemp2+g*Atemp2^3)/(eta*e0*atemp1)
% if dampClass == 1
%     (c*max1+g*max1^3+eta*(e0+w0)*max3)/(eta*(w0+e0)*max3)
%     (c*max3+g*max3^3-eta*(e0+w0)*max1)/(eta*(w0+e0)*max1)
% else
%     (c*max1+g*max1^3+eta*w0*max3)/(eta*w0*max3)
%     (c*max3+g*max3^3-eta*w0*max1)/(eta*w0*max1)
% end




ucNumTot = 45;
ucNumsc = ucNumTot*2;
imagescZ = zeros(size(z,1),4);
for ntemp3 = 1:ucNumsc
    imagescZ(:,ntemp3) = z(1:size(z,1),2*ntemp3-1);
end

figure(110); clf;
hold on;

%imagesc(1-0.5:0.5:ucNumTot-0.5,(1:1:size(z,1))*dt/(2*pi/e0),abs(imagescZ(:,1:2:end)),[0 5])
imagesc(1-0.5:0.5:ucNumTot-0.5,(1:1:size(z,1))*dt/(2*pi/e0),abs(imagescZ(:,1:2:end)),[0 max(max(abs(imagescZ(:,1:2:end))))])
% imagesc(1-0.5:0.5:ucNumTot-0.5,(1:1:size(z,1))*dt/(2*pi/e0),log(abs(imagescZ(:,1:2:end))),[log(3e-1) log(25)])

axis([0 ucNumTot 0 pNum])
xlabel('$n$','Interpreter','latex')
% yticks([-Ac -0.5 0 0.5 Ac])
% yticklabels({'-A_c','-0.5','0','0.5','A_c'})
xticks([15-0.5 30-0.5 45-0.5])
xticklabels([15 30 45])
ylabel('$t/T$','Interpreter','latex')
yticks([0 pNum/2 pNum])
hold off;
colorbar
colormap(hot)
% colorbar('Ticks',[-5,-2,1,4,7],...
%          'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
colorbar('Ticks',[0 Ac 2.5 5])
% legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
% set(gca,'linewidth',2)
% set(gca,'FontSize',19)







sz = 30;
figMag = 1;
lineMag = 10;
mapColor = nan(1001,3);
mapColor2 = nan(1001,4);
SSHucNum = 8;
R1 = 1; 
phi1 = pi/2-pi/6;
r1 = [0 0];
r2 = R1*[cos(phi1) sin(phi1)];
a1 = [2*R1*cos(phi1) 0];
figure(71); clf;
hold on;

bColor = [0 0 255]/256;
rColor = [255 0 0]/256;
kColor = [0 0 0]/256;
%yColor = [255 255 0]/256;
yColor = [255 255 255]/256;

labelMap = 0;
for n = 1:1001
    c1Temp = c2-(c1-c2)/1000*(n-1);
    %colorTemp(1:3) = bColor+(bColor-gColor)/(c1-c2)*(c1Temp-c2)*figMag;
    colorTemp(1:3) = bColor+(bColor-kColor)/(c1-c2)*(c1Temp-c2);
    mapColor(1002-n,1:3) = colorTemp(1:3);
    mapColor2(1002-n,1:3) = colorTemp(1:3);
    mapColor2(1002-n,4) = c1Temp;
end
for n = 1:1001
    c1Temp = c2+(c1-c2)/1000*(n-1);
    %colorTemp(1:3) = bColor-(bColor-rColor)/(c1-c2)*(c1Temp-c2)*figMag;
    colorTemp(1:3) = bColor-(bColor-rColor)/(c1-c2)*(c1Temp-c2);
    mapColor(n+1001,1:3) = colorTemp(1:3);
    mapColor2(n+1001,1:3) = colorTemp(1:3);
    mapColor2(n+1001,4) = c1Temp;
end
for n = 1:1001
    c1Temp = c1+(c1-c2)/1000*(n-1);
    %colorTemp(1:3) = rColor-(rColor-yColor)/(c1-c2)*(c1Temp-c1)*figMag;
    colorTemp(1:3) = rColor-(rColor-yColor)/(c1-c2)*(c1Temp-c1);
    mapColor(n+2002,1:3) = colorTemp(1:3);
    mapColor2(n+2002,1:3) = colorTemp(1:3);
    mapColor2(n+2002,4) = c1Temp;
end

for n = 1:SSHucNum
    X1 = (n-1)*a1(1)+r1(1);
    Y1 = (n-1)*a1(2)+r1(2);
    X2 = (n-1)*a1(1)+r2(1);
    Y2 = (n-1)*a1(2)+r2(2);
    
    if c1n(n) <= c1
        %colorTemp(1:3) = bColor-(bColor-rColor)/(c1-c2)*(c1Temp-c2)*figMag;
        colorTemp1 = rColor-(bColor-rColor)/(c1-c2)*(c1n(n)-c1)*figMag;
    elseif c1n(n) > c1
        %colorTemp(1:3) = rColor-(rColor-yColor)/(c1-c2)*(c1Temp-c1)*figMag;
        colorTemp1 = rColor-(rColor-yColor)/(c1-c2)*(c1n(n)-c1)*figMag;
        %colorTemp1 = rColor+(rColor-gColor)/(c1-c2)*(c1n(n)-c1)*figMag;
    end
    
    if c2n(n) >= c2
        %colorTemp(1:3) = bColor-(bColor-rColor)/(c1-c2)*(c1Temp-c2)*figMag;
        colorTemp2 = bColor-(bColor-rColor)/(c1-c2)*(c2n(n)-c2)*figMag;
    elseif c2n(n) < c2
        %colorTemp(1:3) = bColor+(bColor-gColor)/(c1-c2)*(c1Temp-c2)*figMag;
        colorTemp2 = bColor+(bColor-kColor)/(c1-c2)*(c2n(n)-c2)*figMag;
    end
    

    
    if n ~= SSHucNum 
%         plot([X1 X2],[Y1 Y2],'b-','linewidth',3)
%         plot([X2 X1+a1(1)],[Y2 Y1+a1(2)],'r-','linewidth',3)
        plot([X1 X2],[Y1 Y2],'b-','linewidth',c1n(n)*lineMag,'color',colorTemp1)
        plot([X2 X1+a1(1)],[Y2 Y1+a1(2)],'r-','linewidth',c2n(n)*lineMag,'color',colorTemp2)
    elseif n == SSHucNum
%         plot([X1 X2],[Y1 Y2],'b-','linewidth',3)
%         plot([X2 X2+(X1+a1(1)-X2)/2],[Y2 Y2+(Y1+a1(2)-Y2)/2],'r-','linewidth',3)
%         plot([X2+(X1+a1(1)-X2)/2 X1+a1(1)],[Y2+(Y1+a1(2)-Y2)/2 Y1+a1(2)],'r-.','linewidth',3)
        plot([X1 X2],[Y1 Y2],'b-','linewidth',c1n(n)*lineMag,'color',colorTemp1)
        plot([X2 X2+(X1+a1(1)-X2)/2],[Y2 Y2+(Y1+a1(2)-Y2)/2],'r-','linewidth',c2n(n)*lineMag,'color',colorTemp2)
        plot([X2+(X1+a1(1)-X2)/2 X1+a1(1)],[Y2+(Y1+a1(2)-Y2)/2 Y1+a1(2)],'r-.','linewidth',c2n(n)*lineMag,'color',colorTemp2)
    end
    
    deltay = r2(2)/5;
end
X1 = r1(1);
Y1 = r1(2);
X2 = r2(1);
Y2 = r2(2);
for n = 1:SSHucNum
    X1 = (n-1)*a1(1)+r1(1);
    Y1 = (n-1)*a1(2)+r1(2);
    X2 = (n-1)*a1(1)+r2(1);
    Y2 = (n-1)*a1(2)+r2(2);
    plot(X1,Y1,'b.','MarkerSize',sz*2.5)
    plot(X2,Y2,'r.','MarkerSize',sz*2.5)
end

minc12 = min([c1n;c2n]);
maxc12 = max([c1n;c2n]);
mapColor3 = zeros(2,3);
mapColor4 = zeros(2,4);
label = 0;
for n1 = 1:size(mapColor2,1)
    if mapColor2(n1,4) >= minc12 && mapColor2(n1,4) <= maxc12
        label = label+1;
        mapColor3(label,1:3) = mapColor2(n1,1:3);
        mapColor4(label,1:4) = mapColor2(n1,1:4);
    end
end

colormap(mapColor3)
%colormap(mapColor)
colorbar

[c1Diff,c1MapLabel] = min(abs(mapColor4(:,4)-c1*ones(size(mapColor4,1),1)));
[c2Diff,c2MapLabel] = min(abs(mapColor4(:,4)-c2*ones(size(mapColor4,1),1)));
[minc12Diff,minc12MapLabel] = min(abs(mapColor4(:,4)-minc12*ones(size(mapColor4,1),1)));
[maxc12Diff,maxc12MapLabel] = min(abs(mapColor4(:,4)-maxc12*ones(size(mapColor4,1),1)));

% for n1 = 1:size(mapColor,1)
%     
%     if mapColor2(n1,4)-c2 == 0
%         c2MapLabel = n1;
%     elseif mapColor2(n1,4)-c1 == 0
%         c1MapLabel = n1;
%     elseif mapColor2(n1,4)-minc12 == 0
%         minc12Label = n1;
%     elseif mapColor2(n1,4)-maxc12 == 0
%         maxc12Label = n1;
%     end
% end
colorbar('Ticks',[minc12MapLabel c2MapLabel c1MapLabel maxc12MapLabel]/size(mapColor4,1),'TickLabels',{minc12,c2,c1,maxc12})
%colorbar('Ticks',[minc12 c2 c1 maxc12])

axis equal
axis off
hold off;



toc
