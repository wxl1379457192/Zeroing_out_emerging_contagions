%clear
%clc

data = readtable('E:/zeroCOVID_NPI/Rt_0709/codeV0803/SEIR_1128/dataset.csv',VariableNamingRule='preserve');
NPIefficacy = readtable('E:/zeroCOVID_NPI/Rt_0709/codeV0803/SEIR_1128/NPIefficacy.csv',VariableNamingRule='preserve');
%dataset = testdata;
%load('NPIseffect.mat')
%load('dataset.mat')
%%
%dataset = RtsmoothNPI;
outdir = 'E:/zeroCOVID_NPI/Rt_0709/codeV0803/SEIR_simulation_P3_1129';
%dataset = data(data.citycode == 220100|data.citycode == 310000|data.citycode == 220200,:);
dataset = data(data.citycode == 210800|data.citycode == 220200|data.citycode == 450600| ...
    data.citycode == 210600|data.citycode == 340400| ...
    data.citycode ==140100|data.citycode == 361100|data.citycode == 131000| ...
    data.citycode == 220100|data.citycode == 210100| ...
    data.citycode == 370200|data.citycode == 441900|data.citycode == 110000| ...
    data.citycode == 310000|data.citycode == 320500,:);
%dataset = data(data.citycode == 210800,:);

dataset = dataset(dataset.VG=="omicron",:);
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
%%simulation

for i = 1:size(index,1)
    %for i = 1 : 3
    % dataset for each outbreak
    data = dataset(dataset.city==index(i,1),:);
    data = data(data.original_start == string(index(i,2)),:);
    data = data(data.wave == double(string(index(i,3))),:);
    cases = data.I;
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
    %NPIs_effect_all.m(NPIs_effect_all.parameter=={'Logistics_Management'})];
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
    % + Contact(7) * data.Logistics_Management
    Detect = Detect(1) * data.Mass_screening + Detect(2) * data.Medicine_Management;
    NPIs = [Contact data.Facial_Mask Detect];


    par = ISEIRV_parameter_estimation(data, NPIs);
    R0 = [3,8,13];
    L = [1,4,7];
    all= table();
    for x = 1:length(R0) %??????R0 ?????????3???8, 13
        r = R0(x);
        for y = 1:length(L) %????????????????????????l???
            l= L(y);
            %%%%with different R0
            tspan = length(cases);
            NPI = zeros(tspan,3);
            Controls = zeros(tspan,1);
            s = ISEIRV_R0_model(data, NPI, Controls, par, r, tspan, l);
            %s.Scenario(:) = 'R0 = '+string(r)+'_withoutNPI'+"_Latent"+string(l);
            s.R0(:) = string(r);
            s.("Start.NPI")(:) = '0';
            s.BI(:) = '0';
            s.FI(:) = '0';
            s.SI(:) = '0';
            s.C(:) = '0';
            s.Latent(:) = string(l);

            %%%%with different NPI -timelag ?????????????????????t?????????????????????n???
            %%??????????????????NPI?????????????????????????????????x??????5,10,15(y)??????????????????????????????
            t0 = tic;
            tCount1 = 0;
            all = [all;s];
            for t = 5:5:25
                t1 = tic;
                for n1 = 0:0.1:1
                    %NPIs = [Contact data.Facial_Mask Detect];
                    NPIs = zeros(tspan,3);
                    Con = zeros(tspan,1);
                    %NPIs(:,1) = 0;
                    NPIs(t:end,1) = n1;
                    s = ISEIRV_NPI_model(data, NPIs, Con, par, r, tspan, l);
                    s.R0(:) = string(r);
                    s.("Start.NPI")(:) = string(t);
                    s.BI(:) = string(n1);
                    s.FI(:) = 0;
                    s.SI(:) = 0;
                    s.C(:) = 0;
                    s.Latent(:) = string(l);
                    all = [all;s];
                end
                for n2 = 0:0.1:1
                   % NPIs = [Contact data.Facial_Mask Detect];
                   % NPIs(:,2) = 0;
                    NPIs = zeros(tspan,3);
                    Con = zeros(tspan,1);
                    NPIs(t:end,2) = n2;
                    s = ISEIRV_NPI_model(data, NPIs, Con, par, r, tspan, l);
                    s.R0(:) = string(r);
                    s.("Start.NPI")(:) = string(t);
                    s.BI(:) = 0;
                    s.FI(:) = string(n2);
                    s.SI(:) = 0;
                    s.C(:) = 0;
                    s.Latent(:) = string(l);
                    all = [all;s];
                end
                for n3 = 0:0.1:1
                    %NPIs = [Contact data.Facial_Mask Detect];
                    %NPIs(:,3) = 0;
                    NPIs = zeros(tspan,3);
                    Con = zeros(tspan,1);
                    NPIs(t:end,3) = n3;
                    s = ISEIRV_NPI_model(data, NPIs, Con, par, r, tspan, l);
                    s.R0(:) = string(r);
                    s.("Start.NPI")(:) = string(t);
                    s.BI(:) = 0;
                    s.FI(:) = 0;
                    s.SI(:) = string(n3);
                    s.C(:) = 0;
                    s.Latent(:) = string(l);
                    all = [all;s];
                end
                for c = 0:0.1:1
                    %NPIs = [Contact data.Facial_Mask Detect];
                    %Con = data.Control;
                    NPIs = zeros(tspan,3);
                    Con = zeros(tspan,1);
                    Con(t:end,1) = c;
                    s = ISEIRV_NPI_model(data, NPIs, Con, par, r, tspan, l);
                    s.R0(:) = string(r);
                    s.("Start.NPI")(:) = string(t);
                    s.BI(:) = 0;
                    s.FI(:) = 0;
                    s.SI(:) = 0;
                    s.C(:) = string(c);
                    s.Latent(:) = string(l);
                    all = [all;s];
                end
                tCount1 = tCount1 + toc(t1);
                fprintf(['Time 1 is ', num2str(tCount1), ' seconds.\n'])
                fprintf('R0='+string(r)+"_t="+string(t)+"_Latent"+string(l)+'%% \n')
                %    end
                %end
            end
        end
    end
    writetable(all,[outdir,char('/'+string(data.citycode(1))+"_"+string(i)+'_Scenario.xlsx')],'WriteRowNames',true);
end






