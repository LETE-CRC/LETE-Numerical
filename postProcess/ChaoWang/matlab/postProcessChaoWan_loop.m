%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing script to contruct validation data map following the
% paper of Chao and Wan (2006) IA.
%
% Script structure:
% - Evaporation is calculate in each time step for each band to allow original
% size traceback.
% - Particle data file of converge is read for each time step
%   - traceback is performed considering the time step and the injection
%   stop time, i.e. 1.0s.
% - Arrays are separated according to size bands
% - Variables computation is performed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

dt = 0.25; %time step in which data is recorded
t0 = dt;
tf = 10.0;
t = t0:dt:tf;


%% Group in size bands

%Size bands definition
origBands = [
    1.51	375
    3.03	2300
    6.06	7750
    11.99	12775
    20.2	6925
    28.23	3325
    36.39	1875
    45.07	1100
    62.96	700
    66.42	625
    87.95	475
    113.39	375
    138.58	275
    176.29	250
    227.28	200
    377.78	150
    757.2	100
    ];
% initialize arrays
diam_min = zeros(size(origBands,1),size(t,2));
diam_max = zeros(size(origBands,1),size(t,2));

t_endinj = 1.0;
a=find(t>t_endinj);

% find band limits
for i = 1:size(origBands,1)
    if(i==1)
        diam_min(i,1:(a(1)-1)) = 0.0;
        diam_max(i,1:(a(1)-1)) = (origBands(i+1,1)-origBands(i,1))/2 + origBands(i,1);
    elseif(i==size(origBands,1))
        diam_max(i,1:(a(1)-1)) = (origBands(i,1)-origBands(i-1,1))/2 + origBands(i,1);
        diam_min(i,1:(a(1)-1)) = diam_max(i-1,1:(a(1)-1));
    else
        diam_max(i,1:(a(1)-1)) = (origBands(i+1,1)-origBands(i,1))/2 + origBands(i,1);
        diam_min(i,1:(a(1)-1)) = diam_max(i-1,1:(a(1)-1));
    end
end


%% Read source data
clear data_crude
startRow = 3;
endRow = inf;


for j = 1:size(t,2)
    
    filename = sprintf('/media/sacomano/ac6183b7-5ff6-4f64-995a-6c67c36d7b98/Simulations/Converge/noCOVID19/ValidacaoChao/parcel_data/validationch00000%d.par',j);
    
    data_crude = importParData(filename, startRow, endRow);
    data_crude(:,5) = data_crude(:,5)*2*1e6; % transformation from radius to diameter
    
    % Organize data in size bands and perform variables computations
    % IMPORTANT, datalimits must be adjusted according to the current sampling
    % time, data band msut be specified for each time step
    
    for i = 1:size(origBands,1)
        b              = find(data_crude(:,5)>=diam_min(i,1) & data_crude(:,5)<diam_max(i,1));
        bins(i)        = size(b,1);
        data_band(j,i,1) = origBands(i,1);
        data_band(j,i,2) = bins(i);
        data_band(j,i,3) = mean(data_crude(b,2));    
        data_band(j,i,4) = mean(data_crude(b,1).*data_crude(b,1));
        data_band(j,i,5) = mean(data_crude(b,3).*data_crude(b,3));
        data_band(j,i,6) = mean(data_crude(b,3).*data_crude(b,3) + data_crude(b,3).*data_crude(b,3));
    end
%     M_a = squeeze(data_band(1,:,:))
    clear data_crude
end
%%
%compute band limit sizes for each time step following a predefined
%evaporation model


evap_dt = tf-t_endinj;
u_slip  = 0.0;
t_evap = dt:dt:evap_dt;

diam_out = [];
strSize = size(origBands,1);
tStringSize = size(t_evap,2);
for i = 1:strSize
%     diam_0 = diam_max((strSize+1)-i,1);
    diam_0 = diam_max(i,1);
    % calc evaporation model
    diamString = calcEvap(diam_0,u_slip,t_evap,tStringSize);
    diam_out = [diam_out diamString];
end
    test = 1;






