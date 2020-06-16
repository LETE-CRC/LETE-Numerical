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

dt = 0.10; %time step in which data is recorded
t0 = dt;
tf = 8.0;
t = t0:dt:tf;

evaporation = false;


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
diam_min(1,:) = 0.0;
for i = 1:size(origBands,1)
    if(i==1)
%         diam_min(i,1:(a(1)-1)) = 0.0;
        diam_max(i,:) = (origBands(i+1,1)-origBands(i,1))/2 + origBands(i,1);
    elseif(i==size(origBands,1))
        diam_max(i,:) = (origBands(i,1)-origBands(i-1,1))/2 + origBands(i,1);
        diam_min(i,:) = diam_max(i-1,:);
    else
        diam_max(i,:) = (origBands(i+1,1)-origBands(i,1))/2 + origBands(i,1);
        diam_min(i,:) = diam_max(i-1,:);
    end
end
%% Compute band limit sizes for each time step following a predefined evaporation model
if(evaporation)
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
    clear diamString
end

for k = a(1):size(t,2)    
    diam_max(:,k) = diam_out(k-a(1)+1,:)*1e3;
    for j = 2:size(origBands,1)
        diam_min(j,k) = diam_out(k-a(1)+1,j-1)*1e3;
    end
end
end

%% Read source data
clear data_crude
startRow  = 4;
startRow2 = 2;
endRow    = inf;


for j = 1:size(t,2)
    
%     filename = sprintf('/media/sacomano/ac6183b7-5ff6-4f64-995a-6c67c36d7b98/Simulations/Converge/noCOVID19/ValidacaoChao/parcel_data/validationch00000%d.par',j);
    if j<10
        filename  = sprintf('/media/sacomano/ac6183b7-5ff6-4f64-995a-6c67c36d7b98/Simulations/Fluent/noCOVID19/ValidacaoChao/validacao-uni-grossa-jatoCorrigido/particle/ptcl.mpg000%d',j);
        filename2 = sprintf('/media/sacomano/ac6183b7-5ff6-4f64-995a-6c67c36d7b98/Simulations/Fluent/noCOVID19/ValidacaoChao/validacao-uni-grossa-jatoCorrigido/particle/ptcl.mscl000%d',j);
    elseif j<100
        filename  = sprintf('/media/sacomano/ac6183b7-5ff6-4f64-995a-6c67c36d7b98/Simulations/Fluent/noCOVID19/ValidacaoChao/validacao-uni-grossa-jatoCorrigido/particle/ptcl.mpg00%d',j);
        filename2 = sprintf('/media/sacomano/ac6183b7-5ff6-4f64-995a-6c67c36d7b98/Simulations/Fluent/noCOVID19/ValidacaoChao/validacao-uni-grossa-jatoCorrigido/particle/ptcl.mscl00%d',j);
    end
    
    data_crude = importDataFluentMpg(filename, startRow, endRow);
    diam_crude = importDataFluentMsclNM(filename2, startRow2, endRow);
    B = diam_crude';
    B = B(:)';
    diam_crude = B';
    diam_crude = diam_crude(~isnan(diam_crude));
    data_crude = [data_crude diam_crude];
    data_crude(:,5) = data_crude(:,5)*1e6; % transformation from radius to diameter
    
    % Organize data in size bands and perform variables computations
    % IMPORTANT, datalimits must be adjusted according to the current sampling
    % time, data band must be specified for each time step
    
    for i = 1:size(origBands,1)
%         b              = find(data_crude(:,5)>=diam_min(i,1) & data_crude(:,5)<diam_max(i,1));
        b              = find(data_crude(:,5)>=diam_min(i,j) & data_crude(:,5)<diam_max(i,j));
        bins(i)        = size(b,1);
        data_band(j,i,1) = origBands(i,1);
        data_band(j,i,2) = bins(i);
        data_band(j,i,3) = mean(data_crude(b,3));    
        data_band(j,i,4) = mean(data_crude(b,2).*data_crude(b,2));
        data_band(j,i,5) = mean(data_crude(b,4).*data_crude(b,4));
        data_band(j,i,6) = mean(data_crude(b,2).*data_crude(b,2) + data_crude(b,4).*data_crude(b,4));
    end
%     M_a = squeeze(data_band(1,:,:))
    clear data_crude
end

%% Plot resultswithout consideration of evaporation process

for k = 1:5
    switch k
        case 1
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/012um.dat');
            exp_data1 = importExpData(filename3, 2, endRow);
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/X2_012um.dat');
            exp_X2_data1 = importExpData(filename3, 2, endRow);
        case 2
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/028um.dat');
            exp_data2 = importExpData(filename3, 2, endRow);
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/X2_028um.dat');
            exp_X2_data2 = importExpData(filename3, 2, endRow);
        case 3
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/045um.dat');
            exp_data3 = importExpData(filename3, 2, endRow);
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/X2_045um.dat');
            exp_X2_data3 = importExpData(filename3, 2, endRow);
        case 4
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/087.5um.dat');
            exp_data4 = importExpData(filename3, 2, endRow);
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/X2_087.5um.dat');
            exp_X2_data4 = importExpData(filename3, 2, endRow);
        case 5
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/137.5um.dat');
            exp_data5 = importExpData(filename3, 2, endRow);
            filename3 = sprintf('/home/sacomano/Codes/MATLAB/noCOVID19/expData/X2_137.5um.dat');
            exp_X2_data5 = importExpData(filename3, 2, endRow);
    end
    
end

figure(1)
hold on
plot(t-1,data_band(:,4,3),'Color',[0.1 0.1 .9])
plot(exp_data1(:,1),exp_data1(:,2),'^','MarkerEdgeColor',[0.2 0.2 .9])
plot(t-1,data_band(:,6,3),'Color',[.9 0 0])
plot(exp_data2(:,1),exp_data2(:,2),'s','MarkerEdgeColor',[.9 0 0])
plot(t-1,data_band(:,8,3),'Color',[0 .5 0])
plot(exp_data3(:,1),exp_data3(:,2),'o','MarkerEdgeColor',[0 .5 0])
plot(t-1,data_band(:,11,3),'Color',[.8 .5 0])
plot(exp_data4(:,1),exp_data4(:,2),'x','MarkerEdgeColor',[.8 .5 0])
plot(t-1,data_band(:,13,3),'Color',[0 0 0])
plot(exp_data5(:,1),exp_data5(:,2),'+','MarkerEdgeColor',[0 0 0])
set(gca, 'XScale', 'log')
xlabel(gca,'Time after injection [t(s)]')
ylabel(gca,'Mean vertical position [y(m)]')
hold off
 

figure(2)
hold on
plot(t-1,data_band(:,4,2),'Color',[0.1 0.1 .9])
plot(t-1,data_band(:,6,2),'Color',[.9 0 0])
plot(t-1,data_band(:,8,2),'Color',[0 .5 0])
plot(t-1,data_band(:,11,2),'Color',[.8 .5 0])
plot(t-1,data_band(:,13,2),'Color',[0 0 0])
hold off


figure(3)
hold on
plot(t-1,data_band(:,4,4),'Color',[0.1 0.1 .9])
plot(exp_X2_data1(:,1),exp_X2_data1(:,2),'^','MarkerEdgeColor',[0.2 0.2 .9])
plot(t-1,data_band(:,6,4),'Color',[.9 0 0])
plot(exp_X2_data2(:,1),exp_X2_data2(:,2),'s','MarkerEdgeColor',[.9 0 0])
plot(t-1,data_band(:,8,4),'Color',[0 .5 0])
plot(exp_X2_data3(:,1),exp_X2_data3(:,2),'o','MarkerEdgeColor',[0 .5 0])
plot(t-1,data_band(:,11,4),'Color',[.8 .5 0])
plot(exp_X2_data4(:,1),exp_X2_data4(:,2),'x','MarkerEdgeColor',[.8 .5 0])
plot(t-1,data_band(:,13,4),'Color',[0 0 0])
plot(exp_X2_data5(:,1),exp_X2_data5(:,2),'+','MarkerEdgeColor',[0 0 0])
set(gca, 'YScale', 'log')
xlabel(gca,'Time after injection [t(s)]')
ylabel(gca,'Mean square displacement in X direction [x^2(m^2)]')
hold off

figure(4)
hold on
plot(t-1,data_band(:,4,5),'Color',[0.1 0.1 .9])
plot(exp_X2_data1(:,1),exp_X2_data1(:,2),'^','MarkerEdgeColor',[0.2 0.2 .9])
plot(t-1,data_band(:,6,5),'Color',[.9 0 0])
plot(exp_X2_data2(:,1),exp_X2_data2(:,2),'s','MarkerEdgeColor',[.9 0 0])
plot(t-1,data_band(:,8,5),'Color',[0 .5 0])
plot(exp_X2_data3(:,1),exp_X2_data3(:,2),'o','MarkerEdgeColor',[0 .5 0])
plot(t-1,data_band(:,11,5),'Color',[.8 .5 0])
plot(exp_X2_data4(:,1),exp_X2_data4(:,2),'x','MarkerEdgeColor',[.8 .5 0])
plot(t-1,data_band(:,13,5),'Color',[0 0 0])
plot(exp_X2_data5(:,1),exp_X2_data5(:,2),'+','MarkerEdgeColor',[0 0 0])
set(gca, 'YScale', 'log')
xlabel(gca,'Time after injection [t(s)]')
ylabel(gca,'Mean square displacement in X direction [x^2(m^2)]')
hold off

figure(5)
hold on
plot(t-1,data_band(:,4,6),'Color',[0.1 0.1 .9])
plot(exp_X2_data1(:,1),exp_X2_data1(:,2),'^','MarkerEdgeColor',[0.2 0.2 .9])
plot(t-1,data_band(:,6,6),'Color',[.9 0 0])
plot(exp_X2_data2(:,1),exp_X2_data2(:,2),'s','MarkerEdgeColor',[.9 0 0])
plot(t-1,data_band(:,8,6),'Color',[0 .5 0])
plot(exp_X2_data3(:,1),exp_X2_data3(:,2),'o','MarkerEdgeColor',[0 .5 0])
plot(t-1,data_band(:,11,6),'Color',[.8 .5 0])
plot(exp_X2_data4(:,1),exp_X2_data4(:,2),'x','MarkerEdgeColor',[.8 .5 0])
plot(t-1,data_band(:,13,6),'Color',[0 0 0])
plot(exp_X2_data5(:,1),exp_X2_data5(:,2),'+','MarkerEdgeColor',[0 0 0])
set(gca, 'YScale', 'log')
xlabel(gca,'Time after injection [t(s)]')
ylabel(gca,'Mean square displacement in X direction [x^2(m^2)]')
hold off
%%
% %compute band limit sizes for each time step following a predefined
% %evaporation model
% 
% 
% evap_dt = tf-t_endinj;
% u_slip  = 0.0;
% t_evap = dt:dt:evap_dt;
% 
% diam_out = [];
% strSize = size(origBands,1);
% tStringSize = size(t_evap,2);
% for i = 1:strSize
% %     diam_0 = diam_max((strSize+1)-i,1);
%     diam_0 = diam_max(i,1);
%     % calc evaporation model
%     diamString = calcEvap(diam_0,u_slip,t_evap,tStringSize);
%     diam_out = [diam_out diamString];
% end






