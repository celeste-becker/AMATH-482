%% -----------------TEST 2---------------

%% set up PCA data matrix X - lines (3-10)
clear; close all; clc;

%take same number of frames for each camera (200 frames)
%fill the X matrix with the camera position data for each cam
frameNum = 170;
X = zeros(6,frameNum); 

%% Cam 1 - lines (11 - 31)
load('cam1_2.mat')
numFrames = size(vidFrames1_2,4);
zCam1 = zeros(1,numFrames);
yCam1 = zeros(1,numFrames);

for j = 1:numFrames
    bnwFrame = rgb2gray(vidFrames1_2(200:380,300:400,:,j));   
    [z,y,x] = find(bnwFrame >= 250);

    zCam1(j) = 200 + mean(z); 
    yCam1(j) = 300 + mean(y);
end

%save cam 1 position data to X matrix
X(1,:) = yCam1(1:frameNum);
X(2,:) = zCam1(1:frameNum);

figure(1)
plot(1:frameNum,zCam1(1:frameNum), 'r');

%% Cam 2  - lines (32 - 52)
load('cam2_2.mat')
numFrames2 = size(vidFrames2_2,4);

zCam2 = zeros(numFrames2,1);
yCam2 = zeros(numFrames2,1);
%implay(vidFrames2_1(90:300,260:350,:,:));

for j = 1:numFrames2
    bnwFrame = rgb2gray(vidFrames2_2(:,200:400,:,j));
    [z, y, x] = find(bnwFrame >= 250);
    zCam2(j) = mean(z);
    yCam2(j) = 200 + mean(y);
end

X(3,:) = yCam2(25:frameNum+24);
X(4,:) = zCam2(25:frameNum+24);

hold on
plot(1:frameNum,zCam2(25:frameNum+24), 'b')

%% Cam 3 - lines (53 - 79)
load('cam3_2.mat')
numFrames3 = size(vidFrames3_2,4);

zCam3 = zeros(numFrames3,1);
yCam3 = zeros(1,numFrames3);

for j = 1:numFrames3
    bnwFrame = rgb2gray(vidFrames3_2(:,:,:,j));
    
    [z, y, x] = find(bnwFrame >= 250);
    %because the image is sideways, you need to switch the z and y 
    zCam3(j) = mean(y);
    yCam3(j) = mean(z);
end

%save cam 3 pos data to X matrix
X(5,:) = yCam3(1:frameNum);
X(6,:) = zCam3(1:frameNum);

hold on
plot(1:frameNum,zCam3(1:frameNum),'k')
title('Test 2 - Camera Positions');
legend('Cam 1','Cam 2', 'Cam 3');
ylabel('z position');
xlabel('Number of frames');

%% PCA using SVD - lines (80 - 94)
[m,n]=size(X); % compute data size
mn=mean(X,2); % compute mean for each row
X=X-repmat(mn,1,n); % subtract mean

%svd
[u,s,v]=svd(X'/sqrt(n-1)); % perform the SVD - s = singular values




%energies 
sig=diag(s); %this takes the diagonal values of the singular values
energies = sig.^2/sum(sig.^2);

%% Plot energies - lines (95-103)
figure(5)
% subplot(1,2,1)
plot(energies,'ko','Linewidth',2)
axis([0 25 10^-(18) 1])
ylabel('Energy (log scale)')
xlabel('singular value');
title('TEST 2: Energies using singular values'); 
set(gca,'Fontsize',16,'Xtick',0:5:25)
