%% -----------------TEST 3---------------

%% set up PCA data matrix X - lines (3-10)
clear; close all; clc;

%take same number of frames for each camera (200 frames)
%fill the X matrix with the camera position data for each cam
frameNum = 150;
X = zeros(6,frameNum); 

%% Cam 1 - lines (11 - 30)
load('cam1_3.mat')
numFrames = size(vidFrames1_3,4);
zCam1 = zeros(1,numFrames);
yCam1 = zeros(1,numFrames);
for j = 1:numFrames
    bnwFrame = rgb2gray(vidFrames1_3(230:360,250:350,:,j));   
    [z,y,x] = find(bnwFrame >= 250);
 
    zCam1(j) = 230 + mean(z); 
    yCam1(j) = 250 + mean(y);
end

%save cam 1 position data to X matrix
X(1,:) = yCam1(23:frameNum+22);
X(2,:) = zCam1(23:frameNum+22);

figure(1)
plot(1:frameNum,zCam1(23:frameNum+22), 'r');

%% Cam 2 - lines (31 - 51)
load('cam2_3.mat')
numFrames2 = size(vidFrames2_3,4);

zCam2 = zeros(numFrames2,1);
yCam2 = zeros(numFrames2,1);

for j = 1:numFrames2
    bnwFrame = rgb2gray(vidFrames2_3(150:360,200:400,:,j));
    [z, y, x] = find(bnwFrame >= 250);
    zCam2(j) = 150 + mean(z);
    yCam2(j) = 200 + mean(y);
end

%save cam 2 pos data to X matrix
X(3,:) = yCam2(11:frameNum+10);
X(4,:) = zCam2(11:frameNum+10);

hold on
plot(1:frameNum,zCam2(11:frameNum+10), 'b')

%% Cam 3  - lines (52 - 79)
load('cam3_3.mat')
numFrames3 = size(vidFrames3_3,4);

zCam3 = zeros(numFrames3,1);
yCam3 = zeros(1,numFrames3);

for j = 1:numFrames3
    bnwFrame = rgb2gray(vidFrames3_3(150:350,250:400,:,j));
    
    [z, y, x] = find(bnwFrame >= 235);
    %because the image is sideways, you need to switch the z and y 
    zCam3(j) = 150 + mean(y);
    yCam3(j) = 250 + mean(z);
end

%save cam 3 pos data to X matrix
X(5,:) = yCam3(15:frameNum+14);
X(6,:) = zCam3(15:frameNum+14);

hold on
plot(1:frameNum,zCam3(15:frameNum+14),'k')
title('Test 3 - Camera Positions');
legend('Cam 1','Cam 2', 'Cam 3');
ylabel('z position');
xlabel('Number of frames');

%% PCA using SVD - lines (79 - 91)
[m,n]=size(X); % compute data size
mn=mean(X,2); % compute mean for each row
X=X-repmat(mn,1,n); % subtract mean

%svd
[u,s,v]=svd(X'/sqrt(n-1)); % perform the SVD - s = singular values


%energies 
sig=diag(s); %this takes the diagonal values of the singular values
energies = sig.^2/sum(sig.^2);

%% Plot energies - lines (92-100)
figure(5)
plot(energies,'ko','Linewidth',2)
axis([0 25 10^-(18) 1])
ylabel('Energy (log scale)')
xlabel('singular value');
title('TEST 3: Energies using singular values'); 
set(gca,'Fontsize',16,'Xtick',0:5:25)
