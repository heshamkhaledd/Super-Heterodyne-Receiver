%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Super-Heterodyne Receiver                       %
%                                                                       %
% The purpose of this project is to simulate the basic components of an %
% analog communication system. Specifically, an AM modulator and        %
% a corresponding super-heterodyne receiver.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%% Initializing Channels and Preparing them for Modulation  %%
% Signals Reading %
[BBCArabic_2,BBCArabic_2_samples] = audioread('Short_BBCArabic2.wav');
[FM9090,FM9090_samples] = audioread('Short_FM9090.wav');
[QuranPalestine,QuranPalestine_samples] = audioread('Short_QuranPalestine.wav');
[SkyNewsArabia,SkyNewsArabia_samples] = audioread('Short_SkyNewsArabia.wav');

% Make each signal to a mono-channel one %
BBCArabic_2 = BBCArabic_2(:,1) + BBCArabic_2(:,2);
FM9090 = FM9090(:,1) + FM9090(:,2);
QuranPalestine = QuranPalestine(:,1) + QuranPalestine(:,2);
SkyNewsArabia = SkyNewsArabia(:,1) + SkyNewsArabia(:,2);

% Padding the channels to the maximum size of each one %
Max_Size = max([length(BBCArabic_2),length(FM9090),length(QuranPalestine),length(SkyNewsArabia)]);

BBCArabic_2 = padarray(BBCArabic_2, (Max_Size - length(BBCArabic_2)), 'post');
FM9090 = padarray(FM9090, (Max_Size - length(FM9090)), 'post');
QuranPalestine = padarray(QuranPalestine, (Max_Size - length(QuranPalestine)), 'post');
SkyNewsArabia = padarray(SkyNewsArabia, (Max_Size - length(SkyNewsArabia)), 'post');

% Check if all the sample frequencies are equal for further processing %
if isequal(BBCArabic_2_samples,FM9090_samples,QuranPalestine_samples,SkyNewsArabia_samples)
    Fs = BBCArabic_2_samples;
end
%% Soft Plotting to check %%
subplot(2,2,1);
plot(BBCArabic_2)
title('BBCArabic 2');
xlabel('Time');
subplot(2,2,2);
plot(FM9090);
title('FM9090');
xlabel('Time');
subplot(2,2,3);
plot(QuranPalestine);
title('QuranPalestine');
xlabel('Time');
subplot(2,2,4);
plot(SkyNewsArabia);
title('SkyNewsArabia');
xlabel('Time');
figure;

subplot(2,2,1);
plot(abs(fftshift(fft(BBCArabic_2))));
title('BBCArabic 2');
xlabel('Frequency');
subplot(2,2,2);
plot(abs(fftshift(fft(FM9090))));
title('FM9090');
xlabel('Frequency');
subplot(2,2,3);
plot(abs(fftshift(fft(QuranPalestine))));
title('QuranPalestine');
xlabel('Frequency');
subplot(2,2,4);
plot(abs(fftshift(fft(SkyNewsArabia))));
title('SkyNewsArabia');
xlabel('Frequency');
figure

%% AM Modulation Stage %%

% Generating Carriers with Fc= 100+ n+?F -> where ?F= 50 kHZ %
Carriers_Frequencies = (100000:50000:250000);
n = 1;
% Checking if Nyquist rate critiria is achieved, If not, get the %
% multiplier required to resample the signals %
for Idx = 1:4
    while ((n*Fs/2) < Carriers_Frequencies(Idx) )
        n = n+1;
    end
end
% Resampling the signals with the new sampling rate %
BBCArabic_2 = interp(BBCArabic_2,n);
FM9090 = interp(FM9090,n);
QuranPalestine = interp(QuranPalestine,n);
SkyNewsArabia = interp(SkyNewsArabia,n);

% Update older values to new values after resampling the signals %
Fs = n*Fs;
Max_Size = length(BBCArabic_2);

% Create Carrier's Parameters %
Ts = 1/Fs;
N = 0:1:(Max_Size-1);
Carriers = zeros(4,Max_Size);

for Idx = 1:4
   Carriers(Idx,:) = cos(2*pi*Carriers_Frequencies(Idx)*N*Ts);
end
mod_BBCArabic_2 = BBCArabic_2.*Carriers(1,:)';
mod_FM9090 = FM9090.*Carriers(2,:)';
mod_QuranPalestine = QuranPalestine.*Carriers(3,:)';
mod_SkyNewsArabia = SkyNewsArabia.*Carriers(4,:)';

FDM = mod_BBCArabic_2 + mod_FM9090 + mod_QuranPalestine + mod_SkyNewsArabia;

W_axis = -Fs/2 : (Fs/length(FDM)) : Fs/2 - (Fs/length(FDM));
plot(W_axis,(abs(fftshift(fft(FDM)))));
xlabel('Frequency');
title('FDM Signal');
%% RF Stage %%

