%% -----------------TEST 1---------------

%% set up PCA data matrix X - lines (3-9)
clear; close all; clc;
%take same number of frames for each camera
%fill the X matrix with the camera position data for each cam
frameNum = 150;
X = zeros(6,frameNum); 

%% Cam 1 -  lines (10 - 30) 
load('cam1_1.mat')
numFrames = size(vidFrames1_1,4);
zCam1 = zeros(1,numFrames);
yCam1 = zeros(1,numFrames);
%implay(vidFrames1_1);
for j = 1:numFrames
    bnwFrame = rgb2gray(vidFrames1_1(100:420,300:400,:,j));   
    [z,y,x] = find(bnwFrame >= 250);
 
    zCam1(j) = 100 + mean(z); 
    yCam1(j) = 300 + mean(y);
end

%save cam 1 position data to X matrix
X(1,:) = yCam1(1:frameNum);
X(2,:) = zCam1(1:frameNum);

figure(1)
plot(1:frameNum,zCam1(1:frameNum), 'r');

%% Cam 2 - lines (31 - 51)
load('cam2_1.mat')
numFrames2 = size(vidFrames2_1,4);

zCam2 = zeros(numFrames2,1);
yCam2 = zeros(numFrames2,1);
%implay(vidFrames2_1(90:300,260:350,:,:));

for j = 1:numFrames2
    bnwFrame = rgb2gray(vidFrames2_1(90:300,260:350,:,j));
    [z, y, x] = find(bnwFrame >= 250);
    zCam2(j) = 90 + mean(z);
    yCam2(j) = 260 + mean(y);
end

X(3,:) = yCam2(11:frameNum+10);
X(4,:) = zCam2(11:frameNum+10);

hold on
plot(1:frameNum,zCam2(11:frameNum+10), 'b')

%% Cam 3 - lines (52 - 78)
load('cam3_1.mat')
numFrames3 = size(vidFrames3_1,4);

zCam3 = zeros(numFrames3,1);
yCam3 = zeros(1,numFrames3);

for j = 1:numFrames3
    bnwFrame = rgb2gray(vidFrames3_1(220:290,260:450,:,j));
    
    [z, y, x] = find(bnwFrame >= 250);
    %because the image is sideways, you need to switch the z and y 
    zCam3(j) = 220 + mean(y);
    yCam3(j) = 260 + mean(z);
end

%save cam 3 pos data to X matrix
X(5,:) = yCam3(1:frameNum);
X(6,:) = zCam3(1:frameNum);

hold on
plot(1:frameNum,zCam3(1:frameNum),'k')
title('Test 1 - Camera Positions');
legend('Cam 1','Cam 2', 'Cam 3');
ylabel('z position');
xlabel('Number of frames');

%% PCA using SVD - lines (79 - 91)
[m,n]=size(X); % compute data size
mn=mean(X,2); % compute mean for each row
X=X-repmat(mn,1,n); % subtract mean

%svd
[u,s,v]=svd(X'/sqrt(n-1)); % perform the SVD - s = singular values


%compute energy
sig=diag(s); %this takes the diagonal values of the singular values
energies = sig.^2/sum(sig.^2);

%% Plot energies - lines (92-100)
figure(2)
plot(energies,'ko','Linewidth',2)
axis([0 25 10^-(18) 1])
ylabel('Energy (log scale)')
xlabel('singular value');
title('TEST 1: Energies using singular values'); 
set(gca,'Fontsize',16,'Xtick',0:5:25)
