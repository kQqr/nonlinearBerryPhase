tic

dA = 0.003;
tol = 1e-6;

Alist = dA:dA:0.3;

for n2 = 230:1000
    
    if n2 == 52
        % q = 0;
        filename = strcat('NNM0pp_new');
    elseif n2 == 53
        % q = 1/2;
        filename = strcat('NNM1_2pp_new');
    elseif n2 == 54
        % q = 1/3;
        filename = strcat('NNM1_3pp');
    elseif n2 == 55
        % q = 1/4;
        filename = strcat('NNM1_4pp');
    elseif n2 == 56
        % q = 1/5;
        filename = strcat('NNM1_5pp');
    elseif n2 == 57
        % q = 2/5;
        filename = strcat('NNM2_5pp');
    elseif n2 == 58
        % q = 1/7;
        filename = strcat('NNM1_7pp');
    elseif n2 == 59
        % q = 2/7;
        filename = strcat('NNM2_7pp');
    elseif n2 == 60
        % q = 3/7;
        filename = strcat('NNM3_7pp');
    elseif n2 == 61
        % q = 1/8;
        filename = strcat('NNM1_8pp');
    elseif n2 == 62
        % q = 1/9;
        filename = strcat('NNM1_9pp');
    elseif n2 == 63
        % q = 2/9;
        filename = strcat('NNM2_9pp');
    elseif n2 == 64
        % q = 4/9;
        filename = strcat('NNM4_9pp');
    elseif n2 == 65
        % q = 1/11;
        filename = strcat('NNM1_11pp');
    elseif n2 == 66
        % q = 2/11;
        filename = strcat('NNM2_11pp');
    elseif n2 == 67
        % q = 3/11;
        filename = strcat('NNM3_11pp');
    elseif n2 == 68
        % q = 4/11;
        filename = strcat('NNM4_11pp');
    elseif n2 == 69
        % q = 5/11;
        filename = strcat('NNM5_11pp');
    elseif n2 == 70
        continue;
    elseif n2 == 71
        % q = 2/13;
        filename = strcat('NNM2_13pp');
    elseif n2 == 72
        % q = 3/13;
        filename = strcat('NNM3_13pp');
    elseif n2 == 73
        % q = 4/13;
        filename = strcat('NNM4_13pp');
    elseif n2 == 74
        % q = 5/13;
        filename = strcat('NNM5_13pp');
    elseif n2 == 75
        % q = 6/13;
        filename = strcat('NNM6_13pp');
    elseif n2 == 76
        % q = 1/15;
        filename = strcat('NNM1_15pp');
    elseif n2 == 77
        continue;
    elseif n2 == 78
        % q = 7/15;
        filename = strcat('NNM7_15pp');
    elseif n2 == 79
        % q = abs(1/2-1/3);
        filename = strcat('NNM1_3pm');
    elseif n2 == 80
        % q = abs(1/2-1/5);
        filename = strcat('NNM1_5pm');
    elseif n2 == 81
        % q = abs(1/2-2/5);
        filename = strcat('NNM2_5pm');
    elseif n2 == 82
        % q = abs(1/2-1/7);
        filename = strcat('NNM1_7pm');
    elseif n2 == 83
        % q = abs(1/2-1/8);
        filename = strcat('NNM1_8pm');
    elseif n2 == 84
        % q = abs(1/2-1/9);
        filename = strcat('NNM1_9pm');
    elseif n2 == 85
        % q = abs(1/2-1/16);
        filename = strcat('NNM1_16pm');
    elseif n2 == 86
        % q = abs(1/2-2/7);
        filename = strcat('NNM2_7pm');
    elseif n2 == 87
        % q = abs(1/2-3/7);
        filename = strcat('NNM3_7pm');
    elseif n2 == 88
        % q = abs(1/2-2/9);
        filename = strcat('NNM2_9pm');
    elseif n2 == 89
        % q = abs(1/2-2/13);
        filename = strcat('NNM2_13pm');
    elseif n2 == 90
        % q = abs(1/2-3/11);
        filename = strcat('NNM3_11pm');
    elseif n2 == 91
        % q = abs(1/2-3/20);
        filename = strcat('NNM3_20pm');
    elseif n2 == 92
        % q = abs(1/2-3/29);
        filename = strcat('NNM3_29pm');
    elseif n2 == 93
        % q = abs(1/2-4/29);
        filename = strcat('NNM4_29pm');
    elseif n2 == 94
        % q = abs(1/2-5/24);
        filename = strcat('NNM5_24pm');
    elseif n2 == 95
        % q = abs(1/2-5/29);
        filename = strcat('NNM5_29pm');
    elseif n2 == 96
        % q = abs(1/2-6/29);
        filename = strcat('NNM6_29pm');
    elseif n2 == 97
        % q = abs(1/2-7/16);
        filename = strcat('NNM7_16pm');
    elseif n2 == 98
        % q = abs(1/2-11/29);
        filename = strcat('NNM11_29pm');
    elseif n2 == 99
        continue;
    elseif n2 == 100
        
    elseif n2 == 101
        % q = 4/19;
        filename = strcat('NNM4_19pp');
    elseif n2 == 102
        % q = 5/23;
        filename = strcat('NNM5_23pp');
    elseif n2 == 103
        % q = 7/23;
        filename = strcat('NNM7_23pp');
    elseif n2 == 104
        % q = 8/23;
        filename = strcat('NNM8_23pp');
    elseif n2 == 105
        % q = 9/23;
        filename = strcat('NNM9_23pp');
    elseif n2 == 106
        % q = 10/23;
        filename = strcat('NNM10_23pp');
    elseif n2 == 107
        % q = 2/29;
        filename = strcat('NNM2_29pp');
    elseif n2 == 108
        % q = 6/29;
        filename = strcat('NNM6_29pp');
    elseif n2 == 109
        % q = 8/29;
        filename = strcat('NNM8_29pp');
    elseif n2 == 110
        % q = 9/29;
        filename = strcat('NNM9_29pp');
    elseif n2 == 111
        % q = 10/29;
        filename = strcat('NNM10_29pp');
    elseif n2 == 112
        % q = 9/31;
        filename = strcat('NNM9_31pp');
    elseif n2 == 113
        % q = 10/31;
        filename = strcat('NNM10_31pp');
    elseif n2 == 114
        % q = 11/31;
        filename = strcat('NNM11_31pp');
    elseif n2 == 115
        % q = 12/31;
        filename = strcat('NNM12_31pp');
    elseif n2 == 116
        % q = 13/31;
        filename = strcat('NNM13_31pp');
    elseif n2 == 117
        % q = 14/31;
        filename = strcat('NNM14_31pp');
    elseif n2 == 118
        % q = 14/29;
        filename = strcat('NNM14_29pp');
    elseif n2 == 119
        % q = 15/31;
        filename = strcat('NNM15_31pp');
    elseif n2 == 120
        % q = 13/29;
        filename = strcat('NNM13_29pp');
    elseif n2 == 121
        % q = 17/35;
        filename = strcat('NNM17_35pp');
    elseif n2 == 122
        % q = 25/51;
        filename = strcat('NNM25_51pp');
    elseif n2 == 123
        % q = 37/75;
        filename = strcat('NNM37_75pp');
    elseif n2 == 126
        % q = 9/19;
        filename = strcat('NNM9_19pp');
    elseif n2 == 130
        % q = abs(1/2-1/11);
        filename = strcat('NNM1_11pm');
    elseif n2 == 131
        % q = abs(1/2-1/13);
        filename = strcat('NNM1_13pm');
    elseif n2 == 132
        % q = abs(1/2-1/15);
        filename = strcat('NNM1_15pm');
    elseif n2 == 133
        % q = abs(1/2-1/17);
        filename = strcat('NNM1_17pm');
    elseif n2 == 134
        % q = abs(1/2-1/19);
        filename = strcat('NNM1_19pm');
    elseif n2 == 135
        % q = abs(1/2-1/21);
        filename = strcat('NNM1_21pm');
    elseif n2 == 136
        % q = abs(1/2-1/23);
        filename = strcat('NNM1_23pm');
    elseif n2 == 137
        % q = abs(1/2-1/25);
        filename = strcat('NNM1_25pm');
    elseif n2 == 138
        % q = abs(1/2-1/27);
        filename = strcat('NNM1_27pm');
    elseif n2 == 139
        % q = abs(1/2-1/31);
        filename = strcat('NNM1_31pm');
    elseif n2 == 140
        % q = abs(1/2-1/33);
        filename = strcat('NNM1_33pm');
    elseif n2 == 141
        % q = abs(1/2-1/35);
        filename = strcat('NNM1_35pm');
    elseif n2 == 142
        % q = abs(1/2-1/37);
        filename = strcat('NNM1_37pm');
    elseif n2 == 143
        % q = abs(1/2-1/51);
        filename = strcat('NNM1_51pm');
    elseif n2 == 145
        % q = abs(1/2-4/9);
        filename = strcat('NNM4_9pm');
    elseif n2 ==146
        % q = abs(1/2-5/11);
        filename = strcat('NNM5_11pm');
    elseif n2 == 147
        % q = abs(1/2-9/19);
        filename = strcat('NNM9_19pm');
    elseif n2 == 148
        % q = abs(1/2-10/21);
        filename = strcat('NNM10_21pm');
    elseif n2 == 149
        % q = abs(1/2-11/23);
        filename = strcat('NNM11_23pm');
    elseif n2 == 150
        % q = 1/19;
        filename = strcat('NNM1_19pp');
    elseif n2 == 151
        % q = 2/19;
        filename = strcat('NNM2_19pp');
    elseif n2 == 152
        % q = 3/19;
        filename = strcat('NNM3_19pp');
    elseif n2 == 153
        % q = 5/19;
        filename = strcat('NNM5_19pp');
    elseif n2 == 154
        % q = 6/19;
        filename = strcat('NNM6_19pp');
    elseif n2 == 155
        % q = 1/13;
        filename = strcat('NNM1_13pp');
    elseif n2 == 156
        % q = 1/17;
        filename = strcat('NNM1_17pp');
    elseif n2 == 200
        % q = 65/131;
        filename = strcat('NNM65_131pp');
    elseif n2 == 201
        % q = 45/91;
        filename = strcat('NNM45_91pp');
    elseif n2 == 202
        % q = 30/61;
        filename = strcat('NNM30_61pp');
    elseif n2 == 203
        % q = 27/55;
        filename = strcat('NNM27_55pp');
    elseif n2 == 204
        % q = 2/15;
        filename = strcat('NNM2_15pp');
    elseif n2 == 205
        % q = 4/15;
        filename = strcat('NNM4_15pp');
    elseif n2 == 206
        % q = 2/17;
        filename = strcat('NNM2_17pp');
    elseif n2 == 207
        % q = 3/17;
        filename = strcat('NNM3_17pp');
    elseif n2 == 208
        % q = 4/17;
        filename = strcat('NNM4_17pp');
    elseif n2 == 209
        % q = 5/17;
        filename = strcat('NNM5_17pp');
    elseif n2 == 210
        % q = 6/17;
        filename = strcat('NNM6_17pp');
    elseif n2 == 211
        % q = 7/17;
        filename = strcat('NNM7_17pp');
    elseif n2 == 212
        % q = 8/17;
        filename = strcat('NNM8_17pp');
    elseif n2 == 213
        % q = 7/19;
        filename = strcat('NNM7_19pp');
    elseif n2 == 214
        % q = 8/19;
        filename = strcat('NNM8_19pp');
    elseif n2 == 215
        % q = 1/20;
        filename = strcat('NNM1_20pp');
    elseif n2 == 216
        % q = 3/20;
        filename = strcat('NNM3_20pp');
    elseif n2 == 217
        % q = 1/21;
        filename = strcat('NNM1_21pp');
    elseif n2 == 218
        % q = 2/21;
        filename = strcat('NNM2_21pp');
    elseif n2 == 219
        % q = 4/21;
        filename = strcat('NNM4_21pp');
    elseif n2 == 220
        % q = 5/21;
        filename = strcat('NNM5_21pp');
    elseif n2 == 221
        % q = 8/21;
        filename = strcat('NNM8_21pp');
    elseif n2 == 222
        % q = 10/21;
        filename = strcat('NNM10_21pp');
    elseif n2 == 223
        % q = 1/23;
        filename = strcat('NNM1_23pp');
    elseif n2 == 224
        % q = 2/23;
        filename = strcat('NNM2_23pp');
    elseif n2 == 225
        % q = 3/23;
        filename = strcat('NNM3_23pp');
    elseif n2 == 226
        % q = 4/23;
        filename = strcat('NNM4_23pp');
    elseif n2 == 227
        % q = 6/23;
        filename = strcat('NNM6_23pp');
    elseif n2 == 228
        % q = 11/23;
        filename = strcat('NNM11_23pp');
    elseif n2 == 229
        % q = 1/29;
        filename = strcat('NNM1_29pp');
    elseif n2 == 230
        % q = 3/29;
        filename = strcat('NNM3_29pp');
    elseif n2 == 231
        % q = 4/29;
        filename = strcat('NNM4_29pp');
    elseif n2 == 232
        % q = 5/29;
        filename = strcat('NNM5_29pp');
    elseif n2 == 233
        % q = 7/29;
        filename = strcat('NNM7_29pp');
    elseif n2 == 263
        % q = 1/2-2/11;
        filename = strcat('NNM2_11pm');
    elseif n2 == 264
        % q = 1/2-4/11;
        filename = strcat('NNM4_11pm');
    elseif n2 == 265
        % q = 1/2-6/13;
        filename = strcat('NNM6_13pm');
    elseif n2 == 266
        % q = 1/2-7/15;
        filename = strcat('NNM7_15pm');
    elseif n2 == 267
        % q = 1/2-8/17;
        filename = strcat('NNM8_17pm');
    elseif n2 == 268
        % q = 1/2-1/29;
        filename = strcat('NNM1_29pm');
    elseif n2 == 269
        % q = 1/2-14/29;
        filename = strcat('NNM14_29pm');
    elseif n2 == 279
        % q = 27/56;
        filename = strcat('NNM27_56pp');
    elseif n2 == 280
        % q = 29/60;
        filename = strcat('NNM29_60pp');
    elseif n2 == 281
        % q = 35/72;
        filename = strcat('NNM35_72pp');
    elseif n2 == 282
        % q = 19/39;
        filename = strcat('NNM19_39pp');
    elseif n2 == 283
        % q = 37/76;
        filename = strcat('NNM37_76pp');
    elseif n2 == 284
        % q = 1/2;
        filename = strcat('NNM1_2mp_new');
    else
        continue
    end
    
    filename2 = [filename,'fft'];
    
    load(filename2);
    NNMfftTemp = NNMfft;
    NNMfft = zeros(size(NNMfftTemp,1),size(NNMfftTemp,2));
    
    ASortTemp = nan(size(NNMfftTemp,1)/9,1);
    for n1 = 1:size(NNMfftTemp,1)/9
        ASortTemp(n1) = NNMfftTemp(9*(n1-1)+2,1);
    end
    [ASort,ASortOrder] = sort(ASortTemp);
    for n1 = 1:size(NNMfftTemp,1)/9
        for n985 = 1:9
            NNMfft(9*(n1-1)+n985,:) = NNMfftTemp(9*(ASortOrder(n1)-1)+n985,:);
        end
    end
    
    NNMfftTemp2 = NNMfft;
    NNMfft = zeros(2,size(NNMfftTemp2,2));
    
    labelNNMfft = 0;
    for n986 = 1:size(NNMfftTemp2,1)/9
        if sum(abs(NNMfftTemp2(9*(n986-1)+1,:))) ~= 0
            labelNNMfft = labelNNMfft+1;
            for n987 = 1:9
                NNMfft(9*(labelNNMfft-1)+n987,:) = NNMfftTemp2(9*(n986-1)+n987,:);
            end
        end
    end
    
    nStart = size(NNMfft,1)/9;%size(NNM,1);
    
    NNMlabel = nStart;
    
    q = NNMfft(9*(1-1)+2,7)*2*pi;
    for Ain = Alist
        
        flag = 0;
        
        for n988 = 1:size(NNMfft,1)/9
            Atemp8 = NNMfft(9*(n988-1)+2,1);
            if abs(Atemp8-Ain) < dA
                flag = 1;
                break
            end
        end
        if flag == 0
            
            NNMlabel = NNMlabel+1;
            fftNum = 10;
            startFFTloop = 1;
            randA = 0;
            tolmaxzk = 0.1;
            
            fftOrder = 10;

            tol1 = 1e-6;
            wRecSwitch = 0;
            figSee = 1;
            
            % basic physical parameters of the 1D chain

            e0 = 1.5;%2;
            c1 = 0.25;%0.1;%1
            c2 = NNMfft(2,5)*0.37;%0.35;%2
            
            d = 0;
            d1 = 0.22;%0.22;
            d2 = NNMfft(2,5)*0.02;%0.02;%0.02;
            
            f1 = d1;
            f2 = d2;
            
            Ac = sqrt(-4*(abs(c2)-c1)/3/(abs(d2)-d1));
            
            if abs(q/2/pi) < tol || abs(q/2/pi-1/2) < tol
                ucNum = 2;
            elseif abs(q/2/pi-1/3) < tol || abs(q/2/pi-1/6) < tol
                ucNum = 3;
            else
                ucNum = size(NNMfft,2)/4;
            end
            sw = NNMfft(9*(1-1)+2,4);
            scd = NNMfft(9*(1-1)+2,5);
            tol2 = 1e-2;
            tol3 = 1e-4;
            tol4 = 1e-4;
            
            
            DOF = 4*ucNum;
            timeN = 1000;
            critLoopNum = 10000;
            nu = 1/4;
            
            % initialize the NNM using shooting method
            zk = zeros(timeN+1,DOF);
            gidt = zeros(7,DOF);
            CorrectionList = zeros(DOF,2);
            CorrectionLabel = 0;
            
            zfftTot = nan(timeN,4);
            
            maxzk = NNMfft(9*(1-1)+2,1)*ones(DOF,1); % the peak of the disp
            shiftzk = zeros(timeN+1,DOF); % shift the displacement such that the time starts from the maximum position
            shiftLabel = zeros(DOF,1); % record the shifting label
            shiftT = zeros(DOF,1); % record the shifting time
            IntraLabel = zeros(DOF,1); % record the intracell label shift
            IntraPhase = zeros(DOF,1); % record the intracell phase shift
            InterLabel = zeros(DOF,1); % record the intercell label shift
            InterPhase = zeros(DOF,1); % record the intercell phase shift
            
            tol2 = 1e-2;
            tol3 = 1e-4;
            tol4 = 1e-1;
            tol12 = 1e-6;
            
            c1A = c1+3/4*d1*Ain^2;
            c2A = c2+3/4*d2*Ain^2;
            Delw = sw*sqrt(c1A^2+c2A^2+2*c1A*c2A*cos(q));
            phiq0 = asin(c2A*sin(q)/Delw);
            if abs(cos(phiq0)-(c1A+c2A*cos(q))/Delw) > tol12
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
            % maxzkTemp = Ain;
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
            
            errorRec = nan(critLoopNum+1,1);
            
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
                    zfftTot(timeN*(labelz-1)+1:timeN*labelz,1) = zk(1:timeN,1);
                    zfftTot(timeN*(labelz-1)+1:timeN*labelz,2) = zk(1:timeN,2);
                    zfftTot(timeN*(labelz-1)+1:timeN*labelz,3) = zk(1:timeN,3);
                    zfftTot(timeN*(labelz-1)+1:timeN*labelz,4) = zk(1:timeN,4);
                end
                
                outM = expm(M2)-eye(DOF,DOF);
                
                if loopNum == 1 && figSee == 1
                    [maxzk1,maxzk1Label] = max(zk(1:timeN+1,4*(1-1)+1));
                    cosCompare = zeros(timeN+1,1);
                    for n102 = 1:timeN+1
                        cosCompare(n102) = maxzk1*cos((n102-maxzk1Label)/timeN*2*pi);
                    end
                    
                    maxzk1 = max(zk(:,1));
                end
            end
            
            maxzkEnd = (max(zk(:,1))-min(zk(:,1)))/2;
            
            if abs((maxzkEnd-maxzk(1))/maxzk(1)) < tolmaxzk
                
                n1 = 1;
                [maxzk1,maxzk1Label] = max(zk(1:timeN+1,4*(n1-1)+1));
                cosCompare = zeros(timeN+1,1);
                for n102 = 1:timeN+1
                    cosCompare(n102) = maxzk1*cos((n102-maxzk1Label)/timeN*2*pi);
                end
                if figSee == 1
                    
                    maxzk1 = max(zk(:,1));
                    
                    psiL = zeros(fftOrder,4);
                    
                    for fftLabel = 1:4
                        
                        tol6 = 1e-3;
                        P1x = 0;
                        Q1x = 0;
                        P2x = 0;
                        Q2x = 0;
                        
                        FreqMeasure = 2*pi/dTk/size(zfftTot,1)*(0:size(zfftTot,1)); % the ruler to measure the spectrum: in unit of freqUnit
                        FreqMeasureNorm = FreqMeasure(1:end)/(2*pi/Tk);
                        
                        
                        P2x = 1/size(zfftTot,1)*abs(fft(zfftTot(1:end,fftLabel)));
                        P1x = P2x;
                        P1x = P2x(1:floor(size(P2x,1)/2));
                        P1x(2:end-1) = 2*P1x(2:end-1);
                        
                        [maxP1x,maxP1xLabel] = max(P1x);
                        fftLabelLimit = floor(maxP1xLabel*16);
                        
                        for fftLabel2 = 1:fftOrder
                            [minFreq,minFreqLabel] = min(abs(FreqMeasureNorm-fftLabel2*ones(1,size(FreqMeasure,2))));
                            psiL(fftLabel2,fftLabel) = P1x(minFreqLabel);
                        end
                    end
                end
                
                errorRec2 = errorRec;
                
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
                
                eigM = log(abs(eig(expm(M2))));
                
                psiLRec = zeros(2,4*fftOrder);
                

                for fftLabel = 1:4

                    NNMfft(9*(NNMlabel-1)+1,1:DOF) = zk(1,1:DOF);
                    NNMfft(9*(NNMlabel-1)+2,1) = maxzk(1);
                    NNMfft(9*(NNMlabel-1)+2,2) = 2*pi/Tk;
                    NNMfft(9*(NNMlabel-1)+2,3) = nan; % errorFunc
                    NNMfft(9*(NNMlabel-1)+2,4) = sw;
                    NNMfft(9*(NNMlabel-1)+2,5) = scd;
                    NNMfft(9*(NNMlabel-1)+2,6) = ucNum;
                    NNMfft(9*(NNMlabel-1)+2,7) = q/2/pi;
                    NNMfft(9*(NNMlabel-1)+2,8) = c1;
                    NNMfft(9*(NNMlabel-1)+2,9) = c2;
                    NNMfft(9*(NNMlabel-1)+2,10) = d1;
                    NNMfft(9*(NNMlabel-1)+2,11) = d2;
                    NNMfft(9*(NNMlabel-1)+2,12) = Ac;
                    NNMfft(9*(NNMlabel-1)+2,13) = nan; % stdSum
                    NNMfft(9*(NNMlabel-1)+3,1:DOF) = transpose(InterPhase(1:DOF));
                    NNMfft(9*(NNMlabel-1)+4,1:DOF) = transpose(IntraPhase(1:DOF));
                    NNMfft(9*(NNMlabel-1)+4+fftLabel,1:fftOrder) = transpose(psiL(1:fftOrder,fftLabel)); % Fourier series
                    
                    NNMfft(9*(NNMlabel-1)+2,7) = mean(transpose(InterPhase(1:DOF))); % q
                    NNMfft(9*(NNMlabel-1)+9,1) = (sum(IntraPhase(3:4:DOF))-sum(IntraPhase(1:4:DOF)))/ucNum; %phiq
                    
                end
                
                g2 = 0;
                g2deno = 0;
                g2nume = 0;
                for fftLabel2 = 1:fftOrder
                    g2deno = g2deno+fftLabel2*(psiL(fftLabel2,1)^2+psiL(fftLabel2,2)^2+psiL(fftLabel2,3)^2+psiL(fftLabel2,4)^2);
                    g2nume = g2nume+fftLabel2*(psiL(fftLabel2,1)^2+psiL(fftLabel2,2)^2-psiL(fftLabel2,3)^2-psiL(fftLabel2,4)^2);
                end
                
            end
            
            [B,I] = sort(NNMfft(2:9:end,1));
            NNMfftTemp3 = NNMfft;
            for n989 = 1:size(NNMfftTemp3,1)/9
                for n990 = 1:9
                    NNMfft(9*(n989-1)+n990,:) = NNMfftTemp3(9*(I(n989)-1)+n990,:);
                end
            end
            
            Acheck = NNMfft(2:9:end,1);
            wcheck = NNMfft(2:9:end,2);
            phiqcheck = NNMfft(9:9:end,1);
            qcheck = NNMfft(2:9:end,7);
            
            filename3 = [filename,'fftv2'];
            save(filename3,'NNMfft');
        end
    end
end

toc
