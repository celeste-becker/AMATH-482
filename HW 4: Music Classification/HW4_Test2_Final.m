%% TEST 2
close all; clear all; clc;

%% Get Spectrogram Data

%Chance the Rapper Spectogram data
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest2\Chance\*.mp3');
chanceSpecData = [];

for k = 1:length(songFiles)
        
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Song: #%d = %s\n', k, songFiles(k).name);
    data = spectogramCol(songPath);
    chanceSpecData(:,k) = data;
    
end

%Jcole spectogram data
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest2\JCole\*.mp3');
%Loop Through All Songs in Folder
jcoleSpecData = [];

for k = 1:length(songFiles)
        
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %PRINTS song name
    fprintf('Song: #%d = %s\n', k, songFiles(k).name);
    data = spectogramCol(songPath);
    jcoleSpecData(:,k) = data;
    
end

% Kanye west spectogram data
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest2\Kanye\*.mp3');
%Loop Through All Songs in Folder
kanyeSpecData = [];

for k = 1:length(songFiles)
        
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %PRINTS song name
    fprintf('Song: #%d = %s\n', k, songFiles(k).name);
    data = spectogramCol(songPath);
    kanyeSpecData(:,k) = data;
    
end

%% Training Data
%set the feature number equal to the number of principal components 
%that will be considered
feature = 25;

%Call the trainer function, which does the SVD, and extracts principal
%components, and returns the threshold and the projection of the different
%classes of data onto the w vector
[U,S,V,th1,th2,w,sortChance,sortJCole,sortKanye,musicProjection] = trainer(chanceSpecData, jcoleSpecData, kanyeSpecData, feature);
size(musicProjection)

meanChance = mean(sortChance,2);
meanJCole = mean(sortJCole,2);
meanKanye = mean(sortKanye,2);

%% Test Data - lines (71
%The training data is organized in a folder. Below we go through the
%traing data for each genre and count how many times the program correctly
%classified the song, and then returns the 
songFiles = dir('C:\Users\celes\Documents\UW\AMATH 482\HW\HW 4\DataTest2\TestData\*.mp3');

%JCole test data
songData = [];
counterJCole = 0;
for k = 21:30    
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Test Song: #%d = %s\n', k, songFiles(k).name);
    songData= spectogramCol(songPath);
     
    testSong = U'*songData;
    pval = w'*testSong;
    testClassfications(k) = pval;
    
    %checks if the predicted value for the song is correct for jcole
    if (pval > th1)
        counterJCole = counterJCole+1;
    end
end

 percentCorrectJCole = 100*(counterJCole/10);
 
 
%Kanye West test data
 songData = [];
counterKanye = 0;

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
        counterKanye = counterKanye+1;
    end
end
 percentCorrectKanye = 100*(counterKanye/10);
 
 
%Chance the Rapper Test data
songData = [];
counterChance = 0;
%only  one thing in this for loop so get it to work for one
for k = 1:10    
    %this gets the file path of the song
    songPath = strcat(songFiles(k).folder, "/", songFiles(k).name);

    %prints song name
    fprintf('Test Song: #%d = %s\n', k, songFiles(k).name);
    songData= spectogramCol(songPath);
     
    testSong = U'*songData;
    pval = w'*testSong;
    testClassfications(k) = pval;
    
    if (pval > th2)
        counterChance = counterChance+1;
    end
end

 percentCorrectChance = 100*(counterChance/10);
 
 %calculates the total correct percentage
 totalPercentCorrect = 100*((counterKanye+counterJCole+counterChance)/length(songFiles));
 
%% Print out Percent Accuracy
fprintf('\n');
fprintf('\n');
fprintf('Percent Accuracy for test JCole songs = %d \n', percentCorrectJCole);
fprintf('Percent Accuracy for test Kanye West songs = %d \n', percentCorrectKanye);
fprintf('Percent Accuracy for test Chance songs = %d \n', percentCorrectChance);
fprintf('Total Percent Accuracy = %d \n', totalPercentCorrect);

%% Plot - projections onto w and mean(*)
figure(1)
%chance = blue
plot(sortChance,zeros(),'ob','Linewidth',2);
hold on
m1 = plot(meanChance,0,'*b','Linewidth',10,'DisplayName','jazz mean');

%jcole = red
hold on
plot(sortJCole,1,'or','Linewidth',2);
hold on
m2 = plot(meanJCole,1,'*r','Linewidth',10,'DisplayName','hiphop');

%kanye = black
hold on
plot(sortKanye,ones()*2,'ok','Linewidth',2);
hold on
m3 = plot(meanKanye,2,'*k','Linewidth',10,'DisplayName','metal');

%plot thresholds
hold on
m4 = xline(th1,'g','Linewidth',2);
hold on
xline(th2,'g','Linewidth',2); 

ylim([0 5])
title('Test 2: 1-D Data Projection onto w');
ylabel('Arbitrary number');
xlabel('1-D projection onto w');
legend([m1,m2,m3,m4],'Chance mean','JCole mean','Kanye mean','thresholds')

%% -------- DEFINING FUNCTIONS ------

%Spectogram function. Takes in data and computs the spectrogram, returns
%the data reshaped into one column vector.
function data = spectogramCol(songPath)

        %geting music / spectogram data --> Change song files?
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


%Trainer fucntion - lines (236-311)
%is used on the training data and performs the SVD and LDA
%also computs the first threshold (th1) and second threshold (th2)
function [U,S,V,th1,th2,w,sortgenre1,sortgenre2,sortgenre3, musicProjection] = trainer(genre01,genre02,genre03,feature)

    n1 = size(genre01,2); n2 = size(genre02,2); n3 = size(genre03,2);
    allMusic = [genre01 genre02 genre03];
    
    %SVD
    [U,S,V] = svd(allMusic,'econ');
    
    % PCA projection onto principal components
    musicProjection = S*V'; 
    U = U(:,1:feature);
    Genre1 = musicProjection(1:feature,1:n1);
    Genre2 = musicProjection(1:feature,n1+1:n1+n2);
    Genre3 = musicProjection(1:feature,n2+n3+1:n2+n3+n1);

    m1 = mean(Genre1,2);
    m2 = mean(Genre2,2);
    m3 = mean(Genre3,2);
    mTotal = mean(allMusic(1:feature),2);

    % LDA -- WITHIN CLASS VARIANCES
    Sw = 0; 
    for k=1:n1
        Sw = Sw + (Genre1(:,k)-m1)*(Genre1(:,k)-m1)';
    end
    for k=1:n2
        Sw = Sw + (Genre2(:,k)-m2)*(Genre2(:,k)-m2)';
    end
    for k=1:n3
        Sw = Sw + (Genre3(:,k)-m3)*(Genre3(:,k)-m3)';
    end

    % BETWEEN CLASS VARIANCES  
    SbGenre1 = (m1-mTotal)*(m1-mTotal)';
    SbGenre2 = (m2 - mTotal)*(m2 - mTotal)';
    SbGenre3 = (m3 - mTotal)*(m2 - mTotal)';
    Sb = SbGenre1 + SbGenre2 + SbGenre3; 

    % linear discriminant analysis
    [V2,D] = eig(Sb,Sw); 
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);

    vgenre1 = w'*Genre1; 
    vgenre2 = w'*Genre2;
    vgenre3 = w'*Genre3; 
   
    sortgenre1 = sort(vgenre1);
    sortgenre2 = sort(vgenre2);
    sortgenre3 = sort(vgenre3);
    
    %THRESHOLDS
    %finding threshold 1 - between jazz and hiphop
    t1 = length(sortgenre2);
    t2 = 1;
    while sortgenre2(t1)>sortgenre3(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    %threshold 1
    th1 = (sortgenre2(t1)+sortgenre3(t2))/2;
    
    %finding threshold 2 - between hiphop and metal
    t1 = length(sortgenre3);
    t2 = 1;
    while sortgenre3(t1)>sortgenre1(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    %threshold 2
    th2 = (sortgenre1(t1)+sortgenre3(t2))/2;

end