%% PART 1 Set Up - (lines 1-14)
clear; close all; clc
load handel

amplitude = y'; %amplitude of signal in the time domain
t = (1:length(amplitude))/Fs; %time vector in seconds
S = amplitude; %signal S equal to amplitude

%Defining frequencys and datapoints to eventually make the spectogram 
L=length(t);
n=length(amplitude); % n is the number of data points in signal  
k=(1/L)*[0:(n-1)/2 -(n-1)/2:-1]; %Odd # of frequencys. Frequency in Hz.
ks=fftshift(k); %shifts frequencys to go from negative values to positive

%% Spectrograms for Varying Gabor Widths - (lines 15-46)
a_vec = [0.1 5 50 150]; % a_vec = different filter widths

figure(1)
%outer for loop rotates through the different filter widths 
for jj = 1:length(a_vec)
 
    a = a_vec(jj); %sets new filter width a
    tslide=0:0.1:max(t); %slides the gabor window through the signal in time
    Sgt_spec = zeros(length(tslide),n);
    
    %this computes the spectogram for each different filter width a
    for j=1:length(tslide)
        
        g=exp(-a*(t-tslide(j)).^2); %Gabor Filter
        Sg=g.*S; 
        Sgt=fft(Sg); 
        Sgt_spec(j,:) = fftshift(abs(Sgt)); 
        
    end
    
    %Ploting 
    subplot(2,2,jj)
    pcolor(tslide,ks,Sgt_spec.'), 
    shading interp 
    title(['width a = ',num2str(a)],'Fontsize',16);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    colormap(hot) 
    
end

%% Exploring Oversampling and Under Sampling - (lines 47-92)
%oversampling 
figure(2)
a = 0.001; % this is a  wide filter
tslide=0:0.1:max(t); % slides the time vector
Sgt_spec = zeros(length(tslide),n); %re-defines the spectogram freuqency vector 

for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); 
    Sg=g.*S; 
    Sgt=fft(Sg); 
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 
end

%ploting oversampling spectogram
pcolor(tslide,ks,Sgt_spec.'), 
shading interp 
title('Oversampling Spectogram');
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
set(gca,'Fontsize',16) 
colormap(hot)


%Undersampling
figure(3)
a = 1000; % this is a relly narrow 
tslide=0:1:max(t); % time vector
Sgt_spec = zeros(length(tslide),n); %re-defines the spectogram freuqency vector 

for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); 
    Sg=g.*S; 
    Sgt=fft(Sg); 
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 
end

%ploting under sampling spectogram
pcolor(tslide,ks,Sgt_spec.'), 
shading interp 
title('Under Sampling Spectogram');
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
set(gca,'Fontsize',16) 
colormap(hot)

%% Mexican Hat Wavelet - (lines 93-115)

figure(4)
a = 100; % this is a relatively narrow filter

tslide=0:0.1:max(t); % this makes the time vector "move barely"
Sgt_spec = zeros(length(tslide),n); % this re-defines the spectogram freuqency vector 

for j=1:length(tslide)
    mh = a*(1-(t-tslide(j)).^2).*exp(-((t-tslide(j)).^2)/2); %mexican hat filter
    Sg = mh.*S; 
    Sgt = fft(Sg); 
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 
end

%plot spectogram 
pcolor(tslide,ks,Sgt_spec.'), 
shading interp 
title('Mexican Hat wavelet');
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
set(gca,'Fontsize',16) 
colormap(hot)