%clear
%clc

dataset = readtable('E:/zeroCOVID_NPI/Version0504/dataset/Rt&smooth_NPI.csv',VariableNamingRule='preserve');
NPIefficacy = readtable('E:/zeroCOVID_NPI/Version0504/SEIR/NPIefficacy.csv',VariableNamingRule='preserve');
%dataset = testdata;
%load('NPIseffect.mat')
%load('dataset.mat')
%%
%dataset = RtsmoothNPI;
outdir = 'E:/zeroCOVID_NPI/Version0504/SEIR';
dataset = dataset(dataset.I>=1,:);
%dataset = dataset(dataset.VG == "delta",:);
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

valName  = {'MSE','RE_MSE','r2','VG','citycode','original_start','Predicted_cases','Real_cases'};
valTypes = {'double','double','double','string','string','string','double','double'};
rowcolval = [size(index,1),8];
val = table('Size',rowcolval ,'VariableTypes',valTypes,'VariableNames',valName);

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
    
    % estimation
    switch string(data.VG(1))
        case 'original&alpha'
            [Effect_alpha{i,1},Effect_alpha{i,2},val{i,7},val{i,8},val{i,1},val{i,2},val{i,3}] = ISEIRV(data, NPIs, i, outdir);
            k = mean(Effect_alpha{i,1});
            for j = 1:8
                Par(i,j) = {k(j)};
            end
            val(i,4) = {'Pre-delta'};
            val(i,5) = {data.citycode(1)};
            val(i,6) = {convertTo(data.original_start(1),'yyyymmdd')};
            Par(i,9) = {'Pre-delta'};
            Par(i,10) ={data.citycode(1)};
            Par(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'delta'
            [Effect_delta{i,1},Effect_delta{i,2},val{i,7},val{i,8},val{i,1},val{i,2},val{i,3}] = ISEIRV(data, NPIs, i, outdir);
            k = mean(Effect_delta{i,1});
            for j = 1:8
                Par(i,j) = {k(j)};
            end
            val(i,4) = {'Delta'};
            val(i,5) = {data.citycode(1)};
            val(i,6) = {convertTo(data.original_start(1),'yyyymmdd')};
            %val(i,7) = {data.Label_q(1)};
            Par(i,9) = {'Delta'};
            Par(i,10) ={data.citycode(1)};
            Par(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        case 'omicron'
            [Effect_omicron{i,1},Effect_omicron{i,2},val{i,7},val{i,8},val{i,1},val{i,2},val{i,3}] = ISEIRV(data, NPIs, i, outdir);
            k = mean(Effect_omicron{i,1});
            for j = 1:8
                Par(i,j) = {k(j)};
            end
            val(i,4) = {'Omicron'};
            val(i,5) = {data.citycode(1)};
            val(i,6) = {convertTo(data.original_start(1),'yyyymmdd')};
            %val(i,7) = {data.Label_q(1)};
            Par(i,9) = {'Omicron'};
            Par(i,10) ={data.citycode(1)};
            Par(i,11) ={convertTo(data.original_start(1),'yyyymmdd')};
        otherwise
            fprintf('Error: unexpected variants \n')
    end
    
    

    %
    if ~isempty(Effect_alpha{i,1})
        Results_alpha{i,1} = 1-exp(-Effect_alpha{i,1}(:,1)*mean(NPIs(:,1)));
        Results_alpha{i,2} = 1-exp(-Effect_alpha{i,1}(:,2)*mean(NPIs(:,3)));
        Results(i,:) = Effect_alpha{i,2};

       
    end
    if ~isempty(Effect_delta{i,1})
        Results_delta{i,1} = 1-exp(-Effect_delta{i,1}(:,1)*mean(NPIs(:,1)));
        Results_delta{i,2} = 1-exp(-Effect_delta{i,1}(:,2)*mean(NPIs(:,3)));
        Results(i,:) = Effect_delta{i,2};
    end
    if ~isempty(Effect_omicron{i,1})
        Results_omicron{i,1} = 1-exp(-Effect_omicron{i,1}(:,1)*mean(NPIs(:,1)));
        Results_omicron{i,2} = 1-exp(-Effect_omicron{i,1}(:,2)*mean(NPIs(:,3)));
        Results(i,:) = Effect_omicron{i,2};
    end
    
    fprintf(string(i/size(index,1)*100)+'%% \n')    
end

writetable(Par,[outdir,'/Parameters_SEIR.xlsx'],'WriteRowNames',true)
%writetable(val,[outdir,'/Validation_SEIR.xlsx'],'WriteRowNames',true)
%%
Results_plot = Results(Results(:,1)~=0,:);
Results = Results(Results(:,1)<=1e4,:);
[r, P] = corr(Results);
scatter(log(Results(:,1)),log(Results(:,2)))
hold on
plot(2:max(log(Results(:,1)))+1,2:max(log(Results(:,1)))+1)
xlabel('Estimated total cases number','FontName','Times','FontSize',16)
ylabel('Observed total cases number','FontName','Times','FontSize',16)
annotation('textbox',[0.14,0.82,0.1,0.1],'String','The pairwise linear correlation coefficient = '+string(r(1,2)),'LineStyle','none','FontSize',16,'FontName','Times')
annotation('textbox',[0.34,0.76,0.1,0.1],'String','p = '+string(P(1,2)),'LineStyle','none','FontSize',16,'FontName','Times')


