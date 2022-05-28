tic
% an example of Gaussian excitation in the topological regime

TtoT = 1;

% if TtoT == 1
    fDrive = 7e-2;
    fDriveCoeff = 1;
% else
%     fDrive = 7e-2;
%     fDriveCoeff = 1;
% end

pNum = 5000; % number of periods
LineWidth = 3;
blue = [0 0 256]/256;
red = [256 0 0]/256;
yellow = [255 165 0]/256;
gold1 = [255 215 0]/256;
brown = [139 69 19]/256;
alphaValue = 0.4;

ucNum = 45;
DOF = 4*ucNum;
nu = 1/4;

e0 = 1.5;
dw = 0;
d = 0;
c1 = 0.25;
c2 = 0.37;
% dc = c2-c1;
d1 = 0.22;
d2 = 0.02;
% dd = d2-d1;

disPosi = 7;
disc1 = 0.4;
disc2 = c2;
disd1 = 0.15;
disd2 = d2;




Ac = sqrt(-4*(c2-c1)/3/(d2-d1));

% et = e0+abs(dc);
% eb = e0-abs(dc);



f = zeros(DOF,1);
% in the topo regime, no matter we excite x1, y1, or we excite XN, YN, we always get topo excitation! in the non-topo regime, no matter we excite
% x1, y1, or we excite XN, YN, we always get a small excitation, although also exponential increase as (c2/c1)^{-1} (when c1>c2 non-topo)



timeN = 200000; % 200000


w = e0+dw;
dt = 2*pi/w/(timeN/pNum);
nT = floor(2*pi/w/dt);

tnghw = 3*nT;% Gaussian half width of the driving signal, in units of the time count, nt
tnStart = 15*nT; % starting time of the shaking force of the Gaussian signal

% timeN = floor(600*nT);

eta = 0.000;
tInterval = 10*nT:timeN;


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
            
            if n1 == disPosi
                c1temp = disc1;
                c2temp = disc2;
                d1temp = disd1;
                d2temp = disd2;
            else
                c1temp = c1;
                c2temp = c2;
                d1temp = d1;
                d2temp = d2;
            end
            
            Fa  = e0*yi(n20,4*n1-2)+ d*yi(n20,4*n1-2)^3;
            F1a = c1temp*yi(n20,4*n1+0)+d1temp*yi(n20,4*n1+0)^3;
            if n1 == 1
                F2a = 0;
            else
                F2a = c2temp*yi(n20,4*n1-4)+d2temp*yi(n20,4*n1-4)^3;
            end
            Ga  = +Fa+F1a+F2a-f(4*n1-3)*sin(w*n*dt+phase(4*n1-3))*gFactor-eta*yi(n20,4*n1-3);
            
            Fb  = e0*yi(n20,4*n1-3)+ d*yi(n20,4*n1-3)^3;
            F1b = c1temp*yi(n20,4*n1-1)+d1temp*yi(n20,4*n1-1)^3;
            if n1 == 1
                F2b = 0;
            else
                F2b = c2temp*yi(n20,4*n1-5)+d2temp*yi(n20,4*n1-5)^3;
            end
            Gb  = -Fb-F1b-F2b-f(4*n1-2)*sin(w*n*dt+phase(4*n1-2))*gFactor-eta*yi(n20,4*n1-2);
            
            Fc  = e0*yi(n20,4*n1+0)+ d*yi(n20,4*n1+0)^3;
            F1c = c1temp*yi(n20,4*n1-2)+d1temp*yi(n20,4*n1-2)^3;
            if n1 == ucNum
                F2c = 0;
            else
                F2c = c2temp*yi(n20,4*n1+2)+d2temp*yi(n20,4*n1+2)^3;
            end
            Gc  = +Fc+F1c+F2c-f(4*n1-1)*sin(w*n*dt+phase(4*n1-1))*gFactor-eta*yi(n20,4*n1-1);
            
            Fd  = e0*yi(n20,4*n1-1)+ d*yi(n20,4*n1-1)^3;
            F1d = c1temp*yi(n20,4*n1-3)+d1temp*yi(n20,4*n1-3)^3;
            if n1 == ucNum
                F2d = 0;
            else
                F2d = c2temp*yi(n20,4*n1+1)+d2temp*yi(n20,4*n1+1)^3;
            end
            Gd  = -Fd-F1d-F2d-f(4*n1-0)*sin(w*n*dt+phase(4*n1-0))*gFactor-eta*yi(n20,4*n1-0);

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
%legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('$\Psi_n^{(1)}(t)$','Interpreter','latex')

xticks([0 500/2 500])
yticks([-Ac -0.5 0 0.5 Ac])
yticklabels({'-A_c','-0.5','0','0.5','A_c'})
axis([0 500 -Ac*1.05 Ac*1.05])





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
    plot([n n],[0 ampzMax(2*n-1)],'r-','linewidth',2)
    H2 = plot(n,ampzMax(2*n-0),'bo','markersize',5,'linewidth',2);
    plot([n n],[0 ampzMax(2*n-0)],'b-','linewidth',2)
end
legend([H1 H2],{'${\rm max\,}({\rm Re\,}\Psi_n^{(1)})$','${\rm max\,}({\rm Re\,}\Psi_n^{(2)})$'},'Location','northeast','FontSize',16,'Interpreter','latex');
axis([0 45 0 Ac])
xlabel('Site Labeling','Interpreter','latex')
% ylabel('${\rm max\,}|\Psi_n^{(1)}|$','Interpreter','latex')
xticks([0 15 30 45])
yticks([0 Ac/2 Ac])
yticklabels({'0','A_c/2','A_c'})
hold off;
set(gca,'linewidth',2)
set(gca,'FontSize',19)
% set(gca, 'YScale', 'log')



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
H1 = plot((1:timeN)*dt/(2*pi/e0),fRec(:,1),'b-','linewidth',LineWidth);
% H2 = plot((1:timeN)*dt/(2*pi/e0),fRec(:,2),'r-','linewidth',LineWidth/2);
hold off;
legend([H1],{'${\rm Re\,}S_1(t)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('${\rm Re\,}S_1(t)$','Interpreter','latex')
xticks([0 15 30])
yticks([-fDrive 0 fDrive])
axis([0 30 -fDrive*1.1 fDrive*1.1])





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

bsBot = e0+c1-c2;
bsTop = e0-c1+c2;
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
xticks([2*e0/3 e0 4*e0/3])
yticks([0 0.5 1])
axis([2*e0/3 4*e0/3 0 1])
% set(gca, 'YScale', 'log')










figure(106); clf;
hold on;
H1 = plot((1:size(z,1))*dt/(2*pi/e0),z(:,41),'color',red,'linewidth',LineWidth/3);
H2 = plot((1:size(z,1))*dt/(2*pi/e0),z(:,81),'color',blue,'linewidth',LineWidth/3);
H3 = plot((1:size(z,1))*dt/(2*pi/e0),z(:,121),'color',yellow,'linewidth',LineWidth/3);
hold off;
%legend([H1 H2 H3],{'${\rm Re\,}\Psi_1^{(1)}(t)$','${\rm Re\,}\Psi_2^{(1)}(t)$','${\rm Re\,}\Psi_3^{(1)}(t)$'},'Location','northeast','FontSize',21,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('$\Psi_n^{(1)}(t)$','Interpreter','latex')

xticks([0 500/2 500])
yticks([-0.1 0 0.1])
% yticklabels({'-A_c','-0.5','0','0.5','A_c'})
axis([0 500 -0.2 0.2])






figure(107); clf;
hold on;
H1 = plot([10 10],[10 20],'color',red,'linewidth',LineWidth);
H2 = plot([10 10],[10 20],'color',blue,'linewidth',LineWidth);
H3 = plot([10 10],[10 20],'color',yellow,'linewidth',LineWidth);
hold off;
legend([H1 H2 H3],{'${\rm Re\,}\Psi_{10}^{(1)}$','${\rm Re\,}\Psi_{20}^{(1)}$','${\rm Re\,}\Psi_{30}^{(1)}$'},'Location','northeast','FontSize',21,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$t/T$','Interpreter','latex')
ylabel('$\Psi_n^{(1)}(t)$','Interpreter','latex')
xticks([0 500/2 500])
yticks([-0.2 0 0.2])
% yticklabels({'-A_c','-0.5','0','0.5','A_c'})
axis([0 500 -0.2 0.2])








figure(109); clf;
hold on;

bsBot = e0+c1-c2;
bsTop = e0-c1+c2;
% X = [wNonlinear(1,1); WNonlinear(1,1); WNonlinear(1,1); wNonlinear(1,1)];
X = [bsBot; bsTop; bsTop; bsBot];
Y = [0; 0; 0.012; 0.012];
h = fill(X,Y,gold1);
set(h,'EdgeColor','none')
alpha(alphaValue)

%H1 = plot(FreqMeasure2(1:fftLabelLimit),P1x(1:fftLabelLimit)/max(P1x(1:fftLabelLimit)),'b-','linewidth',LineWidth,'color',brown);
H2 = plot(FreqMeasure(1:fftLabelLimit),ampwz(1:fftLabelLimit,81),'b-','linewidth',3);
% plot([+c1-c2 +c1-c2],[0 3e-4],'k--','linewidth',LineWidth)
% plot([-c1+c2 -c1+c2],[0 3e-4],'k--','linewidth',LineWidth)
hold off;
legend([H2],{'${\rm Re\,}\Psi_{20}^{(1)}(\omega)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
%legend([H1 H2],{'${\rm Re\,}S_1(\omega)$','${\rm Re\,}\Psi_1^{(1)}(\omega)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('$\omega$','Interpreter','latex')
% ylabel('FFT','Interpreter','latex')
xticks([2*e0/3 e0 4*e0/3])
yticks([0 0.005 0.01])
axis([2*e0/3 4*e0/3 0 0.012])
% set(gca, 'YScale', 'log')







[maxA,maxALabel] = max(ampwz(:,1));
% FreqMeasure(maxALabel)
figure(10); clf;
hold on;

% below we look for theoretical fitting
tol = 1e-6;
zeta = (-1+1i*sqrt(3))/2;
topoAmpA = zeros(ucNum,1);
topoAmpA(1) = Aedge;

for n = 1:ucNum-1
    
    if n == disPosi
        c1temp = disc1;
        c2temp = disc2;
        d1temp = disd1;
        d2temp = disd2;
    else
        c1temp = c1;
        c2temp = c2;
        d1temp = d1;
        d2temp = d2;
    end
    
    Apre = topoAmpA(n);
    
    ACube = 1;
    BCube = 0;
    CCube = 4*c2temp/3/d2temp;
    DCube = -d1temp/d2temp*Apre^3-4*c1temp/3/d2temp*Apre;
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
end

H1 = plot(1:ucNum,ampwz(maxALabel,4*(1-1)+1:4:4*(ucNum-1)+1),'r-','linewidth',LineWidth);
% H1 = plot(1:ucNum,ampwzAbs(maxALabel,2*(1-1)+1:2:2*(ucNum-1)+1),'r-','linewidth',LineWidth);
% plot(1:ucNum,ampwz(maxALabel,4*(1-1)+2:4:4*(ucNum-1)+2),'b-')
H2 = plot(1:ucNum,ampwz(maxALabel,4*(1-1)+3:4:4*(ucNum-1)+3),'b-','linewidth',LineWidth);
% H2 = plot(1:ucNum,ampwzAbs(maxALabel,2*(1-1)+2:2:2*(ucNum-1)+2),'b-','linewidth',LineWidth);
H3 = plot(1:ucNum,topoAmpA(:),'k--','linewidth',LineWidth/2);
hold off;
legend([H1 H3 H2],{'${\rm Re\,}\Psi^{(1)}_n(\epsilon_0)$','$\psi^{(1)}_n(\epsilon_0)$','${\rm Re\,}\Psi_n^{(2)}(\epsilon_0)$'},'Location','northeast','FontSize',16,'Interpreter','latex');
% legend([H1 H2],{'$|\Psi^{(1)}_n(\epsilon_0)|$','$|\Psi_n^{(2)}(\epsilon_0)|$'},'Location','northeast','FontSize',16,'Interpreter','latex');
set(gca,'linewidth',2)
set(gca,'FontSize',19)
xlabel('Site labeling','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')
xticks([0 ucNum/3 2*ucNum/3 ucNum])
yticks([1e-6 1e-3 1])
yticklabels({'10^{-6}','10^{-3}','1'})
axis([1 ucNum sqrt(2)*5e-7 1])
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









tInterval0 = 5*nT:timeN;
FreqMeasure = 2*pi/dt/size(tInterval0,2)*(0:(size(tInterval0,2)/2)); % the ruler to measure the spectrum: in unit of freqUnit
tol6 = 1e-3;
P1x = 0;
Q1x = 0;
P2x = 0;
Q2x = 0;

P2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,116)));
P1x = P2x;
P1x = P2x(1:floor(size(tInterval0,2)/2)+1);
P1x(2:end-1) = P1x(2:end-1);

Q2x = 2/size(tInterval0,2)*abs(fft(z(tInterval0,116)));
Q1x = Q2x;
Q1x = Q2x(1:floor(size(tInterval0,2)/2)+1);
Q1x(2:end-1) = Q1x(2:end-1);

[maxP1x,maxP1xLabel] = max(P1x);

figure(12); clf;
hold on;
plot(FreqMeasure(1:min(maxP1xLabel*10,size(FreqMeasure,2))),P1x(1:min(maxP1xLabel*10,size(FreqMeasure,2))),'k-','linewidth',3)
hold off;










tInterval0 = 5*nT:timeN;
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







toc
