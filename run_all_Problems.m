close all;
% clear all;
clc;
%%%%% Problem 1
Delta=0.5;
theta_range=-90:0.1:90;
for M=[2,4,7]
    w = ones(M,1);
    spat_response(w,Delta,theta_range);
    hold on
end
M=7;
w = ones(M,1);
for Delta=[1,2]
    spat_response(w,Delta,theta_range);
end
legend('Delta=0.5,M=2','Delta=0.5,M=4','Delta=0.5,M=7','Delta=1,M=7','Delta=2,M=7');
hold off

%%%% P2
figure;
theta = [5,10];
Delta=0.5;
for M =[7,20]
    for N=[20,100]
        for SNR=[0,20]
            X=gen_data(M,N,Delta,theta,SNR);
            singular_values = svd(X);
            %svd(X)-sort(sqrt(eig((X*X'))),"descend");
            plot(1:M,singular_values,'*');
            title('singular values of X for sep = 5 deg');
            xlabel('singular value index');
            ylabel('singular value');
            hold on
        end
    end
end

legend("theta=[5,10],M=7,N=20,SNR=0","theta=[5,10],M=7,N=20,SNR=20","theta=[5,10],M=7,N=100,SNR=0","theta=[5,10],M=7,N=100,SNR=20","theta=[5,10],M=20,N=20,SNR=0","theta=[5,10],M=20,N=20,SNR=20","theta=[5,10],M=20,N=100,SNR=0","theta=[5,10],M=20,N=100,SNR=20")
hold off

figure;
theta = [5,65];
Delta=0.5;
for M =[7,20]
    for N=[20,100]
        for SNR=[0,20]
            X=gen_data(M,N,Delta,theta,SNR);
            singular_values = svd(X);
            %svd(X)-sort(sqrt(eig((X*X'))),"descend");
            plot(1:M,singular_values,'*');
            title('singular values of X for sep = 60 deg');
            xlabel('singular value index');
            ylabel('singular value');
            hold on
        end
    end
end

legend("theta=[5,65],M=7,N=20,SNR=0","theta=[5,65],M=7,N=20,SNR=20","theta=[5,65],M=7,N=100,SNR=0","theta=[5,65],M=7,N=100,SNR=20","theta=[5,65],M=20,N=20,SNR=0","theta=[5,65],M=20,N=20,SNR=20","theta=[5,65],M=20,N=100,SNR=0","theta=[5,65],M=20,N=100,SNR=20")
hold off

figure;
theta = [5,10];
Delta=0.5;
for M =[7,20]
    for N=[20,100]
        for SNR=[0,20]
            X=gen_data(M,N,Delta,theta,SNR);
            R_hat = 1/N * (X*X');
            singular_values_R = sort(sqrt(eig(N*R_hat)),"descend");
            plot(1:M,singular_values_R,'*');
            title('singular values of X using eig values of R_{hat} for sep = 5 deg');
            xlabel('singular value index');
            ylabel('singular value');
            hold on
        end
    end
end
legend("theta=[5,10],M=7,N=20,SNR=0","theta=[5,10],M=7,N=20,SNR=20","theta=[5,10],M=7,N=100,SNR=0","theta=[5,10],M=7,N=100,SNR=20","theta=[5,10],M=20,N=20,SNR=0","theta=[5,10],M=20,N=20,SNR=20","theta=[5,10],M=20,N=100,SNR=0","theta=[5,10],M=20,N=100,SNR=20")
hold off

figure;
theta = [5,65];
Delta=0.5;
for M =[7,20]
    for N=[20,100]
        for SNR=[0,20]
            X=gen_data(M,N,Delta,theta,SNR);
            R_hat = 1/N * (X*X');
            singular_values_R = sort(sqrt(eig(N*R_hat)),"descend");
            plot(1:M,singular_values_R,'*');
            title('singular values of X using eig values of R_{hat} for sep = 60 deg');
            xlabel('singular value index');
            ylabel('singular value');
            hold on
        end
    end
end
legend("theta=[5,65],M=7,N=20,SNR=0","theta=[5,65],M=7,N=20,SNR=20","theta=[5,65],M=7,N=100,SNR=0","theta=[5,65],M=7,N=100,SNR=20","theta=[5,65],M=20,N=20,SNR=0","theta=[5,65],M=20,N=20,SNR=20","theta=[5,65],M=20,N=100,SNR=0","theta=[5,65],M=20,N=100,SNR=20")
hold off

%%%% Prob 3

% Given: X, A
% To find: S that minimizes ‖ X − AS ‖^2
M=5;N=1000;Delta=0.5;theta=[0,60];SNR=20;
X = gen_data(M,N,Delta,theta,SNR);
% Generate the array response matrix A(theta)
d = length(theta);
A = zeros(M, d);
for i = 1:d
    A(:, i) = gen_a(M, Delta, theta(i));
end

S_est_mf = A' * X;
S_est_zfr = pinv(A) * X;

figure;
plot(S_est_mf(1,:).', 'x')
hold on
plot(S_est_zfr(1,:).', 'o')
% convert the array to a string
array_str = sprintf('[%d,%d] ', theta);

% create a title using the array string
title_str = ['Constellation for theta = ' array_str];
title(title_str);
hold off
legend("Matched filter","Zero Forcing Response")




%%%%%% Prob 4

M=5;N=10;Delta=0.5;theta=[0,60];SNR=20;
X = gen_data(M,N,Delta,theta,SNR);
Rx = 1/N *(X*X');
Rx_inv = Rx^(-1);
P_mvdr = zeros(1,181);
P_c = zeros(1,181);
P_music = zeros(1,181);
%Classical Beamformer
for i = 1:181
    a=gen_a(M,Delta,i-91); % theta from -90 to 90
    P_c(i)=(a'*Rx*a)./(a'*a);
end
% MVDR
for i = 1:181
    a=gen_a(M,Delta,i-91); % theta from -90 to 90
    P_mvdr(i)=(a'*a)/(a'*Rx_inv*a);
end
% MUSIC
var=10^(-SNR/10);
[eigen_vec,eigen_val] = eig(Rx);
Un = eigen_vec(:,1:M-2);
for i = 1:181
    a=gen_a(M,Delta,i-91); % theta from -90 to 90
    P_music(i)=(a'*Un*Un'*a./(1))^(-1);
end
figure
plot(-90:90,10*log10(abs(P_c)));
xlabel('angle[deg]');
ylabel('Power[dB]')

hold on
plot(-90:90,10*log10(abs(P_mvdr)));
plot(-90:90,10*log10(abs(P_music)));
% convert the array to a string
array_str = sprintf('[%d,%d] ', theta);

% create a title using the array string
title_str = ['for theta = ' array_str];
title(title_str);
legend('classical','mvdr','music')
hold off



% (c) function
function X = gen_data(M, N, Delta, theta, SNR)
% Generate a data matrix X = A*theta*S + Nas, where A is the array response matrix,
% S is the source symbols matrix, and Nas is the complex Gaussian noise matrix.

% Generate the array response matrix A(theta)
d = length(theta);
A = zeros(M, d);
for i = 1:d
    A(:, i) = gen_a(M, Delta, theta(i));
end

% Generate the source symbols matrix S
S = (randi([0, 1], d, N)*2-1 - 1j*(randi([0, 1], d, N)*2-1)) / sqrt(2);

% Generate the complex Gaussian noise matrix Nas
sigma = sqrt(10^(-SNR/10)); % calculate noise variance from SNR
Nas = (randn(M, N) + 1j*randn(M, N)) * sigma / sqrt(2);

% Generate the data matrix X = A*theta*S + Nas
X = A * S + Nas;
end


% (b) function
function y = spat_response(w, Delta, theta_range)
    % Define the array response vector
    M = length(w);
    a = exp(1i*2*pi*Delta*(0:M-1)'*sin(deg2rad(theta_range)));
    
    % Compute the spatial response
    y = abs(w'*a);
    
    % Plot the spatial response
    plot(theta_range, y);
    xlabel('angle[deg]');
    ylabel('Spatial response');
    title('Spatial response for fixed w');
end

%(a) function
function a = gen_a(M, Delta, theta)
    % Convert theta from degrees to radians
    theta = deg2rad(theta);
    
    % Construct array response vector
    a = exp(1i*2*pi*Delta*(0:M-1)'*sin(theta));
end
