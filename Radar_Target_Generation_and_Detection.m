clear all
clc;

%% Radar Specifications 
% Frequency of operation = 77GHz
fc = 77e9;
% Max Range = 200m
R_max = 200;
% Range Resolution = 1 m
d_res = 1;
% Max Velocity = 100 m/s
max_velocity = 100;
% speed of light = 3e8
c = 3e8;

%% User Defined Range and Velocity of target
d0 = 80; % Initial position
v0 = -50; % Initial velocity
 
%% FMCW Waveform Generation
% Bandwidth -> lesson 3.1
B_sweep = c/(2*d_res);
% Chirp time -> lesson 3.2
T_chirp = 5.5*2*R_max/c;
% Slop of chirp -> lesson 2.4
slope = B_sweep/T_chirp;   

% The Number of chirp on each sequence
Nd=128;                   
%The number of samples on each chirp. 
Nr=1024;                  
% Total time for samples
t=linspace(0,Nd*T_chirp,Nr*Nd); 

% Create empty vector to store parameters
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

r_t=zeros(1,length(t)); % range covered
td=zeros(1,length(t)); % time delay

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    %For each timestamp update the Range of the Target for constant velocity. 
    r_t(i) = d0 + v0*t(i);
    td(i) = 2*r_t(i)/c;
    %For each time sample, update the transmitted and received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+slope*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + (slope*(t(i)-td(i))^2)/2));
    
    %Mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) .* Rx(i);
    
end

%% RANGE MEASUREMENT
% Perform FFT on the beat signal along the range bins dimension (Nr)
sig_fft = abs(fft(Mix,Nr)./Nr); % FFT + absolute value + normalization
sig_fft = sig_fft(1:(Nr/2)); % Remove symmetric part of FFT signal

% Plot FFT output 
figure('Name','FFT Output Plot'),
plot(sig_fft);
axis ([0 200 0 1]);
xlabel('Measured range');

%% RANGE DOPPLER RESPONSE
% Reshape the vector into Nr*Nd array. 
Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','Range Doppler Map'),
surf(doppler_axis,range_axis,RDM);

%% CFAR implementation
% Define the number of Training Cells in both the dimensions
Tcr = 10;
Tcd = 4;

% Define the number of Guard Cells in both dimensions around the Cell under 
% test (CUT) for accurate estimation
Gcr = 5;
Gcd = 2;

% Offset the threshold by SNR value in dB
offset = 1.4;

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(Nr/2-2*(Tcd+Gcd),Nd-2*(Tcr+Gcr));

% Select the grid that includes the training, guard and test cells
grid_size = (2*Tcr+2*Gcr+1)*(2*Tcd+2*Gcd+1); % lesson 4.5

% Total number of training Cells -> lesson 4.5
num_train_cell = grid_size - (2*Gcr+1)*(2*Gcd+1);

% Create an empty matrix to store detected signal 
CFAR_sig = zeros(size(RDM));

% Design a loop to slides the CUT across range doppler map
for i = 1:(Nr/2-2*(Gcd+Tcd))
    for j = 1:(Nd-2*(Gcr+Tcr))
        
        % Convert the value from logarithmic to linear
        train_cell_patch = db2pow(RDM(i:i+2*(Tcd+Gcd),j:j+2*(Gcr+Tcr)));
        train_cell_patch(Tcd+1:end-Tcd,Tcr+1:end-Tcr) = 0;
        
        % Mean of the sum the signal level with all the training cell
        mean_train_cell = sum(train_cell_patch,'all')/num_train_cell;
        
        % Convert back to logarithmic dB
        noise_level(i,j) = pow2db(mean_train_cell);
        
        % Dynamic threshold of current training patch
        threshold = noise_level(i,j)*offset;
        
        % Filter the noise with dynamic threshold of current patch
        if RDM(i+(Tcd+Gcd),j+(Tcd+Gcr)) > threshold
            CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 1;
        else
            CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 0;
        end
    end
end

% Plot the detected signal after CFAR applied
figure('Name','CA-CFAR Filtered RDM'),
surf(doppler_axis,range_axis,CFAR_sig);
colorbar;