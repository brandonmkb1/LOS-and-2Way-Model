%simulation variables
f = 1E9;
R_start = 100;
R_end = 1000;
N = 10000;
P = 100;
G_TX = 1;
G_RX = 1;
global epsilon_r
epsilon_r = 15;
global h_t
h_t = 10;
global h_r
h_r = 1;
global isTE
isTE = true; 
%calculated constants
global lambda
lambda = 3E8/f;
global d_break
d_break = 4*h_t*h_r/lambda;

%functions for finding loss
function d_d = directPathLength(r)
    global h_t
    global h_r
    d_d = sqrt(r*r + power(h_t-h_r, 2));
end
function R_1 = reflectionPoint(r)
    global h_t
    global h_r
    R_1 = r*h_t/(h_r+h_t);
end
function theta = reflectionAngle(r)
    global h_t
    global h_r
    theta = atan(h_t/reflectionPoint(r));
end 
function [d_1, d_2] = reflectionPathLengths(r)
    global h_t
    global h_r
    d_1 = sqrt(power(reflectionPoint(r),2)+power(h_t, 2));
    d_2 = sqrt(power(r - reflectionPoint(r),2)+power(h_r, 2));
end
function thetal_T = transmissionAngle(r)
    global epsilon_r
    thetal_T = asin(sin(reflectionAngle(r))/sqrt(epsilon_r));
end 
function gamma = reflectionCoefficent(r)
    global epsilon_r
    global isTE
    if isTE
        num = cos(reflectionAngle(r)) - sqrt(epsilon_r)*cos(transmissionAngle(r));
        denom = cos(reflectionAngle(r)) + sqrt(epsilon_r)*cos(transmissionAngle(r));
        gamma = num/denom;
    else 
        num = sqrt(epsilon_r)*cos(reflectionAngle(r)) - cos(transmissionAngle(r));
        denom = sqrt(epsilon_r)*cos(reflectionAngle(r)) + cos(transmissionAngle(r));
        gamma = num/denom;
    end 
end 
function L_p = directPropLoss(r)
    global lambda
    L_p = power(4*pi*r/lambda, 2);
end
function L = loss(r)
    global d_break
    global lambda
    [d_1, d_2] = reflectionPathLengths(r);
    d = directPathLength(r);
    if r > d_break
        bracket = 1 + reflectionCoefficent(r)*d/(d_1 + d_2)*exp(2i*pi/lambda*(d_1+d_2-d));
        L = power(abs(bracket), 2)/directPropLoss(r);
    else
        L = 1/directPropLoss(r);
    end 
end
function dBm = convTodBm(P_r)
    dBm = 10*log10(P_r*1000);
end 

% Section 2 - 100m to 1000m 10 000 users uniform distrobution
r = R_start : (R_end-R_start)/N :R_end; %setting user distribution where we are testing power
direct_1 = []; % getting the values of direct loss
for i = r
    direct_1 = [direct_1, 1/directPropLoss(i)];
end 
twoRay_1 = []; %getting values of 2 ray model
for i = r
    twoRay_1 = [twoRay_1, loss(i)];
end
% generating asked for plots
figure(1)
plot(r,convTodBm(P*G_TX*G_RX*twoRay_1),r,convTodBm(P*G_TX*G_RX*direct_1))
legend("2-ray model", "Direct ray model")
title("Power of received signal with distance")
xlabel("Distance (m)")
ylabel("Power (dBm)")

figure(2)
subplot(2,1,1);
histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','pdf')
ylabel("pdf")
xlabel("Power (dBm)")
title("PDF and CDF of Direct Path Model")
subplot(2,1,2);
histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','cdf')
ylabel("cdf")
xlabel("Power (dBm)")

figure(3)
subplot(2,1,1);
histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','pdf')
ylabel("pdf")
xlabel("Power (dBm)")
title("PDF and CDF of 2-Ray Model")
subplot(2,1,2);
histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','cdf')
ylabel("cdf")
xlabel("Power (dBm)")

%section 3
f = 2e9; % changing frequency, wavelength, d_break to match then plotting
lambda = 3E8/f;
d_break = 4*h_t*h_r/lambda;
direct_1 = [];
for i = r
    direct_1 = [direct_1, 1/directPropLoss(i)];
end 
twoRay_1 = [];
for i = r
    twoRay_1 = [twoRay_1, loss(i)];
end
figure(4)
plot(r,convTodBm(P*G_TX*G_RX*twoRay_1),r,convTodBm(P*G_TX*G_RX*direct_1))
legend("2-ray model", "Direct ray model")
title("Power of received signal with distance for 2GHz")
xlabel("Distance (m)")
ylabel("Power (dBm)")

f = 10e9; % changing f, lambda, d_break again and plotting
lambda = 3E8/f;
d_break = 4*h_t*h_r/lambda;
direct_1 = [];
for i = r
    direct_1 = [direct_1, 1/directPropLoss(i)];
end 
twoRay_1 = [];
for i = r
    twoRay_1 = [twoRay_1, loss(i)];
end
figure(5)
plot(r,convTodBm(P*G_TX*G_RX*twoRay_1),r,convTodBm(P*G_TX*G_RX*direct_1))
legend("2-ray model", "Direct ray model")
title("Power of received signal with distance for 10GHz")
xlabel("Distance (m)")
ylabel("Power (dBm)")

%section 4
f = 1e9; % changing f, lambda, d_break back to the original and plotting
lambda = 3E8/f;
d_break = 4*h_t*h_r/lambda;
isTE = false;
direct_1 = []; % getting the values of direct loss
for i = r
    direct_1 = [direct_1, 1/directPropLoss(i)];
end 
twoRay_1 = []; %getting values of 2 ray model
for i = r
    twoRay_1 = [twoRay_1, loss(i)];
end
% generating asked for plots
figure(6)
plot(r,convTodBm(P*G_TX*G_RX*twoRay_1),r,convTodBm(P*G_TX*G_RX*direct_1))
legend("2-ray model", "Direct ray model")
title("Power of received signal with distance of a TM wave")
xlabel("Distance (m)")
ylabel("Power (dBm)")

figure(7)
subplot(2,1,1);
histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','pdf')
ylabel("pdf")
xlabel("Power (dBm)")
title("PDF and CDF of Direct Path Model of TM wave")
subplot(2,1,2);
histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','cdf')
ylabel("cdf")
xlabel("Power (dBm)")

figure(8)
subplot(2,1,1);
histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','pdf')
ylabel("pdf")
xlabel("Power (dBm)")
title("PDF and CDF of 2-Ray Model of TM wave")
subplot(2,1,2);
histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','cdf')
ylabel("cdf")
xlabel("Power (dBm)")

%section 5
isTE = true;
epsilon_r = 81;

direct_1 = []; % getting the values of direct loss
for i = r
    direct_1 = [direct_1, 1/directPropLoss(i)];
end 
twoRay_1 = []; %getting values of 2 ray model
for i = r
    twoRay_1 = [twoRay_1, loss(i)];
end
% generating asked for plots
figure(9)
plot(r,convTodBm(P*G_TX*G_RX*twoRay_1),r,convTodBm(P*G_TX*G_RX*direct_1))
legend("2-ray model", "Direct ray model")
title("Power of received signal with distance over sea water")
xlabel("Distance (m)")
ylabel("Power (dBm)")

figure(10)
subplot(2,1,1);
histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','pdf')
ylabel("pdf")
xlabel("Power (dBm)")
title("PDF and CDF of Direct Path Model over sea water")
subplot(2,1,2);
histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','cdf')
ylabel("cdf")
xlabel("Power (dBm)")

figure(11)
subplot(2,1,1);
histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','pdf')
ylabel("pdf")
xlabel("Power (dBm)")
title("PDF and CDF of 2-Ray Model over sea water")
subplot(2,1,2);
histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','cdf')
ylabel("cdf")
xlabel("Power (dBm)")

%section 6
R_min = [1 10 100];
count = -1;
for i = R_min
    count = count +1;
    r = i : (R_end-R_start)/N :R_end;
    direct_1 = []; % getting the values of direct loss
    for j = r
        direct_1 = [direct_1, 1/directPropLoss(j)];
    end 
    twoRay_1 = []; %getting values of 2 ray model
    for j = r
       twoRay_1 = [twoRay_1, loss(j)];
    end
    % generating asked for plots
    figure(12+count*9)
    plot(r,convTodBm(P*G_TX*G_RX*twoRay_1),r,convTodBm(P*G_TX*G_RX*direct_1))
    legend("2-ray model", "Direct ray model")
    title("Power of received signal with distance ")
    xlabel("Distance (m)")
    ylabel("Power (dBm)")

    figure(13+count*9)
    subplot(2,1,1);
    histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','pdf')
    ylabel("pdf")
    xlabel("Power (dBm)")
    title("PDF and CDF of Direct Path Model ")
    subplot(2,1,2);
    histogram(convTodBm(P*G_TX*G_RX*direct_1), 'Normalization','cdf')
    ylabel("cdf")
    xlabel("Power (dBm)")

    figure(14+count*9)
    subplot(2,1,1);
    histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','pdf')
    ylabel("pdf")
    xlabel("Power (dBm)")
    title("PDF and CDF of 2-Ray Model ")
    subplot(2,1,2);
    histogram(convTodBm(P*G_TX*G_RX*twoRay_1), 'Normalization','cdf')
    ylabel("cdf")
    xlabel("Power (dBm)")
    pause(1);
end 
