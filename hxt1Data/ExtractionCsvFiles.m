% Visualisation of the data
% change folder to where the data are 
cd('/Users/lucia/Documents/Projects/Hxt_Bay/hxt1Data');

% extract tables of the input and fluorescence data
InputData = readtable('20170302Mock_input.txt','Delimiter','comma')
FluoData = readtable('experiment_20170302Mock_hxt1_data.txt','Delimiter','comma')

figure(1)
plot(InputData{:,'Var1'}*60,InputData{:,'Var2'});
xlabel('Time [min]')
ylabel('Glucose [%w/V]')
ylim([0,1.10])
% I need to have an event based approximation of the input, to speed up
% computation

figure(2)
% Identify the switching times based on the time derivative of the input signal 
deriv = diff(InputData{:,'Var2'});
Indexes = find(deriv,2)+1;
plot(InputData{:,'Var1'}*60,InputData{:,'Var2'}); hold on; 
plot(InputData{Indexes,'Var1'}*60,InputData{Indexes,'Var2'},'o'); hold on; 
xlabel('Time [min]')
ylabel('Glucose [%w/V]')
ylim([0,1.10])

% Alternative approximation
% The second switching time is centered on the mean value of the indexes of
% the non-zero derivatives.
Indexes_all = find(diff(InputData{:,'Var2'}))+1;
Indexes_mean = [Indexes_all(1) round(mean(Indexes_all(2:end)))];
InputApprox_mean = zeros(size(InputData{:,'Var1'}));
InputApprox_mean(Indexes_mean(1):139) = 1;

figure
plot(InputData{:,'Var1'}*60,InputData{:,'Var2'}); hold on; 
plot(InputData{:,'Var1'}*60,InputApprox_mean,'r'); hold on; 
legend('Experimental Input','Approximated Input')
xlabel('Time [min]')
ylabel('Glucose [%w/V]')
ylim([0,1.10])

figure(4)
errorbar(FluoData{:,'Var1'}*60,FluoData{:,'Var2'},FluoData{:,'Var3'});
hold on; 
plot(InputData{:,'Var1'}*60,InputData{:,'Var2'}*15,'r');
hold on; 
plot(InputData{:,'Var1'}*60,InputApprox_mean*15,'k');
xlabel('Time [min]')
ylabel('Fluorescence [AU]')
%% Extraction of csv File for the input (event based representation)
% As a first approximation I would consider the input as a rectangular
% pulse of amplitude 1 and switching times Indexes_mean
GluLev = [0,1,0];
SwitchT = [InputData{1,'Var1'}*60,InputData{Indexes_mean(1),'Var1'}*60,InputData{139,'Var1'}*60];
Glu = GluLev';
SwitchingTimes = SwitchT';

Glu_pre = zeros(length(GluLev),1);
FinalTime = InputData{end,'Var1'}*60.*ones(length(SwitchingTimes),1);

index = linspace(1,length(GluLev),length(GluLev));
rowsi = strread(num2str(index),'%s');
varNames = {'Switchingtimes','FinalTime','Glupre','Glu'};
Data_table= table(SwitchingTimes,FinalTime,Glu_pre,Glu,'RowNames',rowsi,'VariableNames',varNames);
cd('/Users/lucia/Documents/Projects/Hxt_Bay/ExperimentalDatacsv')

writetable(Data_table,strcat('20170302Hxt1','_Events_Inputs.csv'));   

%% Extraction of csv File for the observable

time_gfp_min = FluoData{:,'Var1'}*60;
gfp_mean = FluoData{:,'Var2'};
gfp_std = FluoData{:,'Var3'};

index = linspace(1,length(time_gfp_min),length(time_gfp_min));
rowsi = strread(num2str(index),'%s');
varNames = {'timeGFP','GFPmean','GFPstd'};
Data_rep= table(time_gfp_min,gfp_mean,gfp_std,'RowNames',rowsi,'VariableNames',varNames);
cd('/Users/lucia/Documents/Projects/Hxt_Bay/ExperimentalDatacsv')
writetable(Data_rep,strcat('20170302Hxt1','_Observables.csv'));    


