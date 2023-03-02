%clear
%clc

data = readtable('E:/zeroCOVID_NPI/Rt_0709/codeV0803/SEIR_1128/dataset.csv',VariableNamingRule='preserve');
NPIefficacy = readtable('E:/zeroCOVID_NPI/Rt_0709/codeV0803/SEIR_1128/NPIefficacy.csv',VariableNamingRule='preserve');
%dataset = testdata;
%load('NPIseffect.mat')
%load('dataset.mat')
%%
%dataset = RtsmoothNPI;
outdir = 'E:/zeroCOVID_NPI/Rt_0709/codeV0803/SEIR_simulation_P4_1129';
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
    for x = 1:length(R0) %设定R0 分别为3，6，9，12,15
        r = R0(x);
        for y = 1:length(L) %设定潜隐期分别为l天
            l= L(y);
            %%%%NPI combined: 假设Facial mask一直存在0.5：实施CC：0.6，0.7;
            %%%%SD:0,0.25,0.5,0.75, 分别从第7、14天开始实施

            tspan = 100;
            NPI = zeros(tspan,3);
            Controls = zeros(tspan,1);

            %%%%with different NPI -timelag 考虑在疫情发生t日后升级管控至n级
            %%每个情境中，NPI都在疫情出现下降后的第x天以5,10,15(y)天为等差数列逐渐放开
            t0 = tic;
            tCount1 = 0;
            for t = 7:7:14
                t1 = tic;
                for cn = 0.6:0.1:0.8 %CT intensity
                    for n1 = 0:0.25:1 % SD强度
                        %NPIs = [Contact data.Facial_Mask Detect];
                        NPIs = zeros(tspan,3);
                        Con = zeros(tspan,1);
                        %NPIs(:,1) = 0;
                        NPIs(t:end,1) = n1;
                        NPIs(t:end,2) = 0.5;
                        Con(t:end,1) = cn;
                        s = ISEIRV_NPIcombined_model(data, NPIs, Con, par, r, tspan, l);
                        s.R0(:) = string(r);
                        s.("Start.NPI")(:) = string(t);
                        s.BI(:) = string(n1);
                        s.FI(:) = 0.5;
                        s.SI(:) = 0;
                        s.C(:) = string(cn);
                        s.Latent(:) = string(l);
                        all = [all;s];
                    end
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






