%% Stalink satellite downlink design
% 
%  CREATED BY:             Person Code:           Matriculation Number:
%    Niccolò Bardazzi          10800456                         963039
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  SATELLITE COMMUNICATION  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%         AND         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  POSITIONING SYSTEMS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear  
clc
addpath(genpath(pwd))

% constants
k = 1.38e-23;        % Boltzmann's constant
c = 3e8;             % [m/s] speed of light
Tph = 313.15;        % [K] maximum achievable physical temperature
T0 = 290;            % [K] reference temperature

% INPUT variables
h = 505.4e3;         % [m] altitude
f = 12.7e9;          % [Hz] frequency
R0 = 1e9;            % 1 Gbps, data-rate 
Drx = 0.48;          % [m], dish diameter
P_av = 0.99;         % = 99%
lat = 45.5;          % [°] Milan latitude
lon = 9.2;           % [°] Milan longitude
el_angle = 90;       % [°] elevation angle
B = 500e6;           % 250MHz, bandwidth

% Components' parameters
A_wg = 0.8;          % linear attenuation of the waveguide
NFdB = 1.9;          % [dB] LNA noise figure

% Common engineering values:
e = 1;               % [°] pointing error
HPBW = 5.5;          % [°] Half Power Beamwidth
etarx = 0.567;       % RX antenna's efficiency
etatx = 0.6;         % TX antenna's efficiency

% Parameters
lambda = c/f;

%% Attenuation losses
% free-space losses
Lfs = 10*log10((c/(f*4*pi*h))^2);   

% misaligment losses
Lmis = -12*(e/HPBW)^2;

% Outage probability
P_out = 1-P_av;            % Time percentage of excess for the total  
cfg = p618Config('TotalAnnualExceedance', P_out, 'AntennaDiameter', Drx, ...
    'AntennaEfficiency', etarx, 'Frequency', f, 'Latitude', lat, ...
    'Longitude', lon, 'PolarizationTiltAngle', 90,'ElevationAngle',90); 
[pl,~,T_sky] =  p618PropagationLosses(cfg);
Latm = -pl.At;

% Total Losses
LtotdB = Latm+Lmis+Lfs;
Ltot = 10^(LtotdB/10);

% Uncomment for outage probability vs attenuation to be considered
% n = 100;
% attenuation = zeros(n,1);
% P_out_vec = linspace(0.001,0.999,n);
% for i = 1:n
%     cfg.TotalAnnualExceedance = 1-P_out_vec(i);  % Time percentage of excess for the total 
%     [pl,xpd,tsky] =  p618PropagationLosses(cfg);
%     attenuation(i) = pl.At;
% end
% figure()
% plot((1-P_out_vec)*100,attenuation)
% title('System reliability', 'Interpreter','latex')
% xlabel('Outage probability [\%]', 'Interpreter','latex')
% ylabel('$L_{atm}$ [dB]','Interpreter','latex')
% axis tight
% grid on

%% Reciever side
Grx = (pi*Drx/lambda)^2*etarx;  
GrxdB = 10*log10(Grx);

% System Noise
% Waveguide
T_wg = (1-A_wg)*Tph;

% LNA
NF = 10^(NFdB/10);
T_LNA = (NF-1)*T0;

% Total
T_sys = T_sky+T_wg+T_LNA;

%% Modulation 
% up to 64 QAM scheme chosen by Starlink
BER = 1.e-6;      % BER selection
EbNo_margin = 2;  % dB
EbNo = (0:0.01:25)';

ber4 = berawgn(EbNo,'qam',4);
ber8 = berawgn(EbNo,'qam',8);
ber16 = berawgn(EbNo,'qam',16);
ber32 = berawgn(EbNo,'qam',32);
ber64 = berawgn(EbNo,'qam',64);

figure()
semilogy(EbNo,[ber4 ber8 ber16 ber32 ber64],'LineWidth',2), hold on
xlabel('Eb/No (dB)', 'Interpreter','latex')
ylabel('BER', 'Interpreter','latex')
legend('4-QAM','8-QAM','16-QAM','32-QAM','64-QAM')
title("Theoretical Bit Error Rate",'Interpreter','latex')
grid on
ylim([1.e-7 1])

[~, n4] = min(abs(BER-ber4));
EbN04 = EbNo(n4)+EbNo_margin;
[~, n8] = min(abs(BER-ber8));
EbN08 = EbNo(n8)+EbNo_margin;
[~, n16] = min(abs(BER-ber16));
EbN016 = EbNo(n16)+EbNo_margin;
[~, n32] = min(abs(BER-ber32));
EbN032 = EbNo(n32)+EbNo_margin;
[~, n64] = min(abs(BER-ber64));
EbN064 = EbNo(n64)+EbNo_margin;
scatter(EbNo(n4), ber4(n4), 'r', 'DisplayName', 'Nominal condition','MarkerFaceColor','flat')

%% SNR minimum for all the modulation schemes
R4 = R0/2;        M4 = log2(4);
R8 = R0/3;        M8 = log2(8);
R16 = R0/4;       M16 = log2(16);
R32 = R0/5;       M32 = log2(32);
R64 = R0/6;       M64 = log2(64);

EsN04 = EbN04+10*log10(M4);
EsN08 = EbN08+10*log10(M8);
EsN016 = EbN016+10*log10(M16);
EsN032 = EbN032+10*log10(M32);
EsN064 = EbN064+10*log10(M64);

SNR4 = EsN04+10*log10(R4)-10*log10(B);
SNR8 = EsN08+10*log10(R8)-10*log10(B);
SNR16 = EsN016+10*log10(R16)-10*log10(B);
SNR32 = EsN032+10*log10(R32)-10*log10(B);
SNR64 = EsN064+10*log10(R64)-10*log10(B);

%% Link Budget & Transmitter sizing
% Minimum of SNR --> 4-QAM --> Nominal working point
SNRmin = SNR4;
Prx = 10^(SNRmin/10)*k*T_sys*B;    % linear
PrxdB = 10*log10(Prx);
EIRP = PrxdB-LtotdB-GrxdB;

CNdB = EsN04+10*log10(f)-10*log10(B);

Power = 20;
Gtx = 10^(EIRP/10)/Power;
GtxdB = 10*log10(Gtx);
Dtx = sqrt(Gtx/etatx)*lambda/pi;
fprintf('The transmitter diameter is: %.3f m \n', Dtx)





