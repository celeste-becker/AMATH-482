%% TEST 1
close all; clear all; clc;

%% Get Spectrogram Data - lines (4-60)

%Jazz Data - Duke Ellington
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest1\DukeEllington\*.mp3');
jazzSpecData = [];

%loop through song files in folder
for k = 1:length(songFiles)
        
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Song: #%d = %s\n', k, songFiles(k).name);
    
    %stores spectrogram data for each song in a column 
    jazzSpecData(:,k) = spectogramCol(songPath);
    
end


%Hip Hop Data - J. Cole
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest1\JCole\*.mp3');
hiphopSpecData = [];

%loop through song files in folder
for k = 1:length(songFiles)
        
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Song: #%d = %s\n', k, songFiles(k).name);
    
    %stores spectrogram data for each song in a column 
    hiphopSpecData(:,k) = spectogramCol(songPath);
    
end


%Metal Data - Metallica
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest1\Metalica\*.mp3');
metalSpecData = [];
%loop through song files in folder
for k = 1:length(songFiles)
        
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Song: #%d = %s\n', k, songFiles(k).name);
    
    %stores spectrogram data for each song in a column
    metalSpecData(:,k) = spectogramCol(songPath);
    
end

%% Training Data - lines (61-74)
%set the feature number equal to the number of principal components 
%that will be considered
feature = 24;

%Call the trainer function, which does the SVD, and extracts principal
%components, and returns the threshold and the projection of the different
%classes of data onto the w vector
[U,S,V,th1,th2,w,sortJazz,sortHiphop,sortMetal] = trainer(jazzSpecData, hiphopSpecData, metalSpecData, feature);

meanJazz = mean(sortJazz,2);
meanHiphop = mean(sortHiphop,2);
meanMetal = mean(sortMetal,2);

%% Test Data lines(75-149)
%The training data is organized in a folder. Below we go through the
%traing data for each genre and count how many times the program correctly
%classified the song, and then returns the 
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest1\TestData\*.mp3');

%Metal Test Data
songData = [];
counterMetal = 0;
%Loop Through All Songs in Folder
for k = 1:10    
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Test Song: #%d = %s\n', k, songFiles(k).name);
    songData= spectogramCol(songPath);
     
    testSong = U'*songData;
    pval = w'*testSong;
    testClassfications(k) = pval;
    
    if (pval > th1)
        counterMetal = counterMetal+1;
    end
end
 percentCorrectMetal = 100*(counterMetal/10);
 
 
%Jazz Test Data
songData = [];
counterJazz = 0;
for k = 11:20    
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Test Song: #%d = %s\n', k, songFiles(k).name);
    songData= spectogramCol(songPath);
     
    testSong = U'*songData;
    pval = w'*testSong;
    testClassfications(k) = pval;
    
    if ((pval > th1) && (pval < th2))
        counterJazz = counterJazz+1;
    end
end
 percentCorrectJazz = 100*(counterJazz/10);
 
 
%Hip Hop Test Data
songData = [];
counterHiphop = 0;
%only  one thing in this for loop so get it to work for one
for k = 21:30    
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Test Song: #%d = %s\n', k, songFiles(k).name);
    songData= spectogramCol(songPath);
     
    testSong = U'*songData;
    pval = w'*testSong;
    testClassfications(k) = pval;
    
    if (pval > th2)
        counterHiphop = counterHiphop+1;
    end
end

 percentCorrectHiphop = 100*(counterHiphop/10);
 totalPercentCorrect = 100*((counterJazz+counterMetal+counterHiphop)/length(songFiles));
 
 %% Print out Percent Accuracy - lines (150-157)
fprintf('\n');
fprintf('\n');
fprintf('Percent Accuracy for test Metal songs = %d \n', percentCorrectMetal);
fprintf('Percent Accuracy for test Jazz songs = %d \n', percentCorrectJazz);
fprintf('Percent Accuracy for test Hip Hop songs = %d \n', percentCorrectHiphop);
fprintf('Total Percent Accuracy = %d \n', totalPercentCorrect);

%% Plot - projections onto w and mean(*)- lines(158 - 188)
figure(1)
%jazz = blue
plot(sortJazz,zeros(),'ob','Linewidth',2);
hold on
m1 = plot(meanJazz,0,'*b','Linewidth',10,'DisplayName','jazz mean');

%hiphop = red
hold on
plot(sortHiphop,1,'or','Linewidth',2);
hold on
m2 = plot(meanHiphop,1,'*r','Linewidth',10,'DisplayName','hiphop');

%metal = black
hold on
plot(sortMetal,ones()*2,'ok','Linewidth',2);
hold on
m3 = plot(meanMetal,2,'*k','Linewidth',10,'DisplayName','metal');

%plot thresholds
hold on
m4 = xline(th1,'g','Linewidth',2);
hold on
xline(th2,'g','Linewidth',2); 

ylim([0 5])
title('Test 1: 1-D Data Projection onto w');
ylabel('Arbitrary number');
xlabel('1-D projection onto w');
legend([m1,m2,m3,m4],'Jazz mean','Hiphop mean','Metal mean','thresholds')

%% -------- DEFINING FUNCTIONS ------ lines ( 189-314)

%Spectogram function - lines (191-234)
%Takes in data and computs the spectrogram, returns
%the data reshaped into one column vector.
function data = spectogramCol(songPath)

        %geting music / spectogram data 
        [y,Fs] = audioread(songPath);
        
         %Lower sampling rate by a factor of 10
         Fs = Fs/10;
         song = y(1:Fs:5,:);        
         S = song'; %Signal S equal to the song data amplitude  
        L=5;
        time=length(song)/Fs; % total recorded time in seconds
        %time=length(L)/Fs; % total recorded time in seconds
        %t = (1:L)/Fs; % vector of time in seconds
        
        n=Fs*L; % n is the number of data points in signal  
        t = linspace(0,L,n);
        
        
        k=(1/L)*[0:n/2-1 -n/2:-1]; %EVEN # of frequencys in Hz
        ks=fftshift(k);

        %Take the Spectogram of the song
        tslide=0:0.1:L; %Descrete time vector
        Sgt_spec = zeros(length(tslide),n);
        
        %gabor filter width
        a = 150;
    
        for jj=1:length(tslide)
            %applying gabor filter
            g=exp(-a*(t-tslide(jj)).^2); 
            Sg=g.*S; 
            Sgt=fft(Sg);%fourier transforming to the freuqency domain

            %storing all the filtered signal at that time in the spectogram vector
            Sgt_spec(jj,:) = fftshift(abs(Sgt)); 

        end 

        data = reshape(Sgt_spec, (length(tslide)*length(Sgt_spec)), 1);
end



%Trainer fucntion 
%is used on the training data and performs the SVD and LDA
%also computs the first threshold (th1) and second threshold (th2)
function [U,S,V,th1,th2,w,sortJazz,sortHiphop,sortMetal] = trainer(jazz0,hiphop0,metal0,feature)


    nj = size(jazz0,2); nh = size(hiphop0,2); nm = size(metal0,2);
    allMusic = [jazz0 hiphop0 metal0];
    
    %SVD
    [U,S,V] = svd(allMusic,'econ');
    
    % PCA projection onto principal components
    musicProjection = S*V'; 
    U = U(:,1:feature);
    Jazz = musicProjection(1:feature,1:nj);
    Hiphop = musicProjection(1:feature,nj+1:nj+nh);
    Metal = musicProjection(1:feature,nh+nm+1:nh+nm+nj);

    mj = mean(Jazz,2);
    mh = mean(Hiphop,2);
    mm = mean(Metal,2);
    mTotal = mean(allMusic(1:feature),2);

    % LDA -- WITHIN CLASS VARIANCES
    Sw = 0; 
    for k=1:nj
        Sw = Sw + (Jazz(:,k)-mj)*(Jazz(:,k)-mj)';
    end
    for k=1:nh
        Sw = Sw + (Hiphop(:,k)-mh)*(Hiphop(:,k)-mh)';
    end
    for k=1:nm
        Sw = Sw + (Metal(:,k)-mm)*(Metal(:,k)-mm)';
    end

    % BETWEEN CLASS VARIANCES  
    SbJazz = (mj-mTotal)*(mj-mTotal)';
    SbHiphop = (mh - mTotal)*(mh - mTotal)';
    SbMetal = (mm - mTotal)*(mh - mTotal)';
    Sb = SbJazz + SbHiphop + SbMetal; 


    %linear discriminant analysis
    [V2,D] = eig(Sb,Sw); 
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);

    vjazz = w'*Jazz; 
    vhiphop = w'*Hiphop;
    vmetal = w'*Metal; 

    sortJazz = sort(vjazz);
    sortHiphop = sort(vhiphop);
    sortMetal = sort(vmetal);
      
    %THRESHOLDS
    %finding threshold 1 - between metal and jazz
    t1 = length(sortMetal);
    t2 = 1;
    while sortMetal(t1)>sortJazz(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    %threshold 1
    th1 = (sortMetal(t1)+sortJazz(t2))/2;
    
    %finding threshold 2 - between hiphop and metal
    t1 = length(sortJazz);
    t2 = 1;
    while sortJazz(t1)>sortHiphop(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    %threshold 2
    th2 = (sortJazz(t1)+sortHiphop(t2))/2;   

end