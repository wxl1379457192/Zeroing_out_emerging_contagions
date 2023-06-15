%clear
%clc
dataset = readtable('E:/zeroCOVID_NPI/Version0504/dataset/Rt&smooth_NPI.csv',VariableNamingRule='preserve');
NPIefficacy = readtable('E:/zeroCOVID_NPI/Version0504/SEIR/NPIefficacy.csv',VariableNamingRule='preserve');
%dataset = testdata;
%load('NPIseffect.mat')
%load('dataset.mat')
%%
%dataset = RtsmoothNPI;
outdir = 'E:/zeroCOVID_NPI/Version0504/SEIR_SE';
dataset = dataset(dataset.I>=1,:);
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

varNames  = {'b0','b1','r0','r1','r2','r3','E0','I0','variant','ID','Os'};
varTypes = {'double','double','double','double','double','double','double','double','string','string','string'};
rowcol = [size(index,1),11];
Par = table('Size',rowcol ,'VariableTypes',varTypes,'VariableNames',varNames);
Par1 = table('Size',rowcol ,'VariableTypes',varTypes,'VariableNames',varNames);
Par2 = table('Size',rowcol ,'VariableTypes',varTypes,'VariableNames',varNames);
Par3 = table('Size',rowcol ,'VariableTypes',varTypes,'VariableNames',varNames);
for i = 1 : size(index,1)
%for i = 1 : 10
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
       % + Contact(7) * data.Logistics_Management;
    Detect = Detect(1) * data.Mass_screening + Detect(2) * data.Medicine_Management;
    NPIs = [Contact data.Facial_Mask Detect];
    % SE for DIS
    switch string(data.VG(1))
        case 'original&alpha'
            [Effect_alpha{i,1},Effect_alpha{i,2}] = ISEIRV_DIS_model(data, NPIs);
            k = mean(Effect_alpha{i,1});
            for j = 1:8
                Par(i,j) = {k(j)};
            end
            Par(i,9) = {string(data.VG(1))};
            Par(i,10) ={data.citycode(1)};
            Par(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'delta'
            [Effect_delta{i,1},Effect_delta{i,2}] =  ISEIRV_DIS_model(data, NPIs);
            k = mean(Effect_delta{i,1});
            for j = 1:8
                Par(i,j) = {k(j)};
            end
            Par(i,9) = {string(data.VG(1))};
            Par(i,10) ={data.citycode(1)};
            Par(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'omicron'
            [Effect_omicron{i,1},Effect_omicron{i,2}] = ISEIRV_DIS_model(data, NPIs);
            k = mean(Effect_omicron{i,1});
            for j = 1:8
                Par(i,j) = {k(j)};
            end
            Par(i,9) = {string(data.VG(1))};
            Par(i,10) ={data.citycode(1)};
            Par(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        otherwise
            fprintf('Error: unexpected variants \n')
    end

    % SE for IS1
    switch string(data.VG(1))
        case 'original&alpha'
            [Effect_alpha{i,1},Effect_alpha{i,2}] = ISEIRV_IS1_model(data, NPIs);
            k = mean(Effect_alpha{i,1});
            for j = 1:8
                Par1(i,j) = {k(j)};
            end
            Par1(i,9) = {string(data.VG(1))};
            Par1(i,10) ={data.citycode(1)};
            Par1(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'delta'
            [Effect_delta{i,1},Effect_delta{i,2}] =  ISEIRV_IS1_model(data, NPIs);
            k = mean(Effect_delta{i,1});
            for j = 1:8
                Par1(i,j) = {k(j)};
            end
            Par1(i,9) = {string(data.VG(1))};
            Par1(i,10) ={data.citycode(1)};
            Par1(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'omicron'
            [Effect_omicron{i,1},Effect_omicron{i,2}] =  ISEIRV_IS1_model(data, NPIs);
            k = mean(Effect_omicron{i,1});
            for j = 1:8
                Par1(i,j) = {k(j)};
            end
            Par1(i,9) = {string(data.VG(1))};
            Par1(i,10) ={data.citycode(1)};
            Par1(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        otherwise
            fprintf('Error: unexpected variants \n')
    end
   
    % SE for IS2
    switch string(data.VG(1))
        case 'original&alpha'
            [Effect_alpha{i,1},Effect_alpha{i,2}] = ISEIRV_IS2_model(data, NPIs);
            k = mean(Effect_alpha{i,1});
            for j = 1:8
                Par2(i,j) = {k(j)};
            end
            Par2(i,9) = {string(data.VG(1))};
            Par2(i,10) ={data.citycode(1)};
            Par2(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'delta'
            [Effect_delta{i,1},Effect_delta{i,2}] = ISEIRV_IS2_model(data, NPIs);
            k = mean(Effect_delta{i,1});
            for j = 1:8
                Par2(i,j) = {k(j)};
            end
            Par2(i,9) = {string(data.VG(1))};
            Par2(i,10) ={data.citycode(1)};
            Par2(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'omicron'
            [Effect_omicron{i,1},Effect_omicron{i,2}] = ISEIRV_IS2_model(data, NPIs);
            k = mean(Effect_omicron{i,1});
            for j = 1:8
                Par2(i,j) = {k(j)};
            end
            Par2(i,9) = {string(data.VG(1))};
            Par2(i,10) ={data.citycode(1)};
            Par2(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        otherwise
            fprintf('Error: unexpected variants \n')
    end

    % SE for IS3
    switch string(data.VG(1))
        case 'original&alpha'
            [Effect_alpha{i,1},Effect_alpha{i,2}] = ISEIRV_IS3_model(data, NPIs);
            k = mean(Effect_alpha{i,1});
            for j = 1:8
                Par3(i,j) = {k(j)};
            end
            Par3(i,9) = {string(data.VG(1))};
            Par3(i,10) ={data.citycode(1)};
            Par3(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'delta'
            [Effect_delta{i,1},Effect_delta{i,2}] = ISEIRV_IS3_model(data, NPIs);
            k = mean(Effect_delta{i,1});
            for j = 1:8
                Par3(i,j) = {k(j)};
            end
            Par3(i,9) = {string(data.VG(1))};
            Par3(i,10) ={data.citycode(1)};
            Par3(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'omicron'
            [Effect_omicron{i,1},Effect_omicron{i,2}] = ISEIRV_IS3_model(data, NPIs);
            k = mean(Effect_omicron{i,1});
            for j = 1:8
                Par3(i,j) = {k(j)};
            end
            Par3(i,9) = {string(data.VG(1))};
            Par3(i,10) ={data.citycode(1)};
            Par3(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        otherwise
            fprintf('Error: unexpected variants \n')
    end
   
   
    fprintf(string(i/size(index,1)*100)+'%% \n')    
end

writetable(Par,[outdir,'/DIS_Parameters_SEIR.xlsx'],'WriteRowNames',true)
writetable(Par1,[outdir,'/IS1_Parameters_SEIR.xlsx'],'WriteRowNames',true)
writetable(Par2,[outdir,'/IS2_Parameters_SEIR.xlsx'],'WriteRowNames',true)
writetable(Par3,[outdir,'/IS3_Parameters_SEIR.xlsx'],'WriteRowNames',true)

