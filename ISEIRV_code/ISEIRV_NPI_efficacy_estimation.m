%clear
%clc

dataset = readtable('E:/zeroCOVID_NPI/Version0504/dataset/Rt&smooth_NPI.csv',VariableNamingRule='preserve');
NPIefficacy = readtable('E:/zeroCOVID_NPI/Version0504/SEIR/NPIefficacy.csv',VariableNamingRule='preserve');
%dataset = testdata;
%load('NPIseffect.mat')
%load('dataset.mat')
%%
%dataset = RtsmoothNPI;
outdir = 'E:/zeroCOVID_NPI/Version0504/SEIR_P2';
%dataset = data(data.citycode == 220100|data.citycode == 310000|data.citycode == 340200,:);
index = [dataset.city...
    categorical(dataset.original_start)...
    categorical(dataset.wave)];
index = unique(index,'rows');

Effect_alpha = cell(size(index,1),2);
Effect_delta = cell(size(index,1),2);
Effect_omicron = cell(size(index,1),2);

Results_alpha = Effect_alpha;
Results_delta = Effect_delta;
Results_omicron = Effect_omicron;
Results = zeros(length(Effect_alpha),2);

Par = zeros(length(Effect_alpha),11);
Par = array2table(Par);
Par.Properties.VariableNames(1:11) = {'b0','b1','r0','r1','r2','r3','E0','I0','variant','ID','Os'};

val = zeros(size(index,1),6);
val = array2table(val);
val.Properties.VariableNames(1:6) = {'MSE','RE_MSE','r2','VG','citycode','original_start'};

for i = 1 : size(index,1)
%for i = 1 : 3
    % dataset for each outbreak
    data = dataset(dataset.city==index(i,1),:);
    data = data(data.original_start == string(index(i,2)),:);
    data = data(data.wave == double(string(index(i,3))),:);
    
    if height(data)<=5
        continue
    end

    % aggregate NPIs
    NPIs_effect_all = NPIefficacy(strcmp(NPIefficacy.variant,data.VG(1))==1,:);
    Contact = [NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'Lockdown'})==1)...
        NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'Business_Premises_Closure'})==1)...
        NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'Public_Transportation_Closure'})==1)...
        NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'Gathering_restriction'})==1)...
        NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'Workplace_Closure'})==1)...
        NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'School_Closure'})==1)];
    Contact = Contact/sum(Contact);
    Detect = [NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'Mass_screening'})==1)...
        NPIs_effect_all.m(strcmp(NPIs_effect_all.parameter,{'Medicine_Management'})==1)];
    Detect = Detect/sum(Detect);
    Contact = Contact(1) * data.Lockdown...
        + Contact(2) * data.Business_Premises_Closure...
        + Contact(3) * data.Public_Transportation_Closure...
        + Contact(4) * data.Gathering_restriction...
        + Contact(5) * data.Workplace_Closure...
        + Contact(6) * data.School_Closure;
    Detect = Detect(1) * data.Mass_screening + Detect(2) * data.Medicine_Management;
    NPIs = [Contact data.Facial_Mask Detect];
    
    par = ISEIRV_parameter_estimation(data, NPIs);
    
    
    %%%%real world
    tspan = length(data.I);
    NPI = smoothdata([NPIs(1,:); NPIs],'gaussian',28);
    Controls = data.Control;
    Controls = smoothdata(Controls,'gaussian',14);
    real = ISEIRV_realworld_model(data, NPI,Controls, par, tspan);
    real(:,6) = {'realworld'};
    
    %%%%without any NPI
    tspan = length(data.I);
    NPI = zeros(tspan,3);
    Controls = zeros(tspan,1);
    base = ISEIRV_realworld_model(data, NPI,Controls, par, tspan);
    base(:,6) = {'Baseline'};

    %%%without Contact tracing%%%%
    tspan = length(data.I);
    NPI = smoothdata([NPIs(1,:); NPIs],'gaussian',28);
    Controls = zeros(tspan,1);
    sm1 = ISEIRV_realworld_model(data, NPI,Controls, par, tspan);
    sm1(:,6) = {'withoutCT'};

    %%%without social distance%%%%
    tspan = length(data.I);
    NPI = smoothdata([NPIs(1,:); NPIs],'gaussian',28);
    NPI(:,1) = 0;
    %%NPI(:,3) = 0;
    Controls = data.Control;
    Controls = smoothdata(Controls,'gaussian',14);
    sm2 = ISEIRV_realworld_model(data, NPI,Controls, par, tspan);
    sm2(:,6) = {'withoutSD'};

    %%%without facial mask%%%%
    tspan = length(data.I);
    NPI = smoothdata([NPIs(1,:); NPIs],'gaussian',28);
    NPI(:,2) = 0;
    Controls = data.Control;
    Controls = smoothdata(Controls,'gaussian',14);
    sm3 = ISEIRV_realworld_model(data, NPI, Controls, par, tspan);
    sm3(:,6) = {'withoutFM'};
 
    %%%without PCR%%%%
    tspan = length(data.I);
    NPI = smoothdata([NPIs(1,:); NPIs],'gaussian',28);
    NPI(:,3) = 0;
    Controls = data.Control;
    Controls = smoothdata(Controls,'gaussian',14);
    sm4 = ISEIRV_realworld_model(data, NPI, Controls, par, tspan);
    sm4(:,6) = {'withoutPCR'};

    all = [real; base; sm1; sm2; sm3; sm4];
    
    %%%%with different R0
    %tspan = 1500;
    %NPI = zeros(tspan,3);
    %NPI(1:length(NPIs),:) = NPIs(1:length(NPIs),:); 
    %NPI(length(NPIs)+1:end,1) = NPIs(length(NPIs),1); %Contact
    %NPI(length(NPIs)+1:end,2) = NPIs(length(NPIs),2); %Fical mask
    %NPI(length(NPIs)+1:end,3) = NPIs(length(NPIs),3); %detect: PCR screening
    %Controls = zeros(tspan,1);
    %Controls(1:length(Control),1) = Control(1:length(Control),1); 
    %Controls(length(Control)+1:end,1) = Control(length(Control),1); 

    %R0list = {3,6,9,12};
    %for r = 1:length(R0list)
    %  s = SEIR_Simulation_R0(data, NPI, Controls, par, R0list(r), tspan);
    %  s(:,6) = {string('R0 = ',R0list(r))};
    %  all = [all;s];
    %end
    writetable(all,[outdir,char('/Scenario_'+string(i)+'.xlsx')],'WriteRowNames',true)
    %[Effect_omicron{i,1}, Effect_omicron{i,2}] = SEIR_Simulation_S1(data, NPIs, outdir, 1.25,1.5,0.75,0.5);
    %[Effect_omicron{i,1}, Effect_omicron{i,2}] = SEIR_Simulation_S2(data, NPIs, outdir, 1.0, 3.0, 5.0, 7.0);
    %[Effect_omicron{i,1}, Effect_omicron{i,2}] = SEIR_Simulation_S3(data, NPIs, outdir, -0.25, -0.5, 0.25, 0.5);
    %[Effect_omicron{i,1}, Effect_omicron{i,2}] = SEIR_Simulation_withoutNPI(data, NPIs, outdir);
    fprintf(string(i/size(index,1)*100)+'%% \n')    
end

%index2 = [dataset.ID];
%index2 = unique(index2,'rows');
%for i = 1:size(index2,1)
%    data = dataset(dataset.ID==index2(i,1),:);
    %[Effect_omicron{i,1}, Effect_omicron{i,2}] = SEIR_Simulation_S4(data,NPIefficacy, outdir);
    %[Effect_omicron{i,1}, Effect_omicron{i,2}] = SEIR_Simulation_S4v2(data,NPIefficacy, outdir);
    %[Effect_omicron{i,1}, Effect_omicron{i,2}] = SEIR_Simulation_S4v3(data,NPIefficacy, outdir);
%end

