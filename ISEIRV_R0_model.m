function [sim] = ISEIRV_R0_model(data, NPI, Controls, par, R0, tspan, latent)

% parameter initialization
pop = data.pop(1)*1e4;        % population
vac = pop*data.("Practical.vaccination");
vac(2:end) = vac(2:end) - vac(1:end-1);
%%

b0 = R0/5;
% sampling
%b0 = makedist('normal','mu',b0,'sigma',1);
l = 1/latent; % 1/ latent period
%% output_simulation
samples = 4000;

b1 = par(2);
r0 = par(3);
r1 = par(4);
r2 = par(5);
r3 = par(6);
E0 = par(7);
I0 = par(8);

I_new = zeros(tspan,samples);
I_det = zeros(tspan,samples);
S = zeros(tspan,1);
E = zeros(tspan,1);
I = zeros(tspan,1);
r = zeros(tspan,samples);
b = zeros(tspan,samples);
variable_num = 8;
Effect = zeros(samples,variable_num);
for ss = 1 : samples
    k0 = [b0 b1.random r0.random r1.random r2.random r3.random E0.random I0.random];
    %
    Effect(ss,:) = k0;
    %
    for t = 1 : tspan
        b(t,ss) = k0(1)*exp(-k0(2)*NPI(t,1))*(1-0.25*NPI(t,2));
        if (k0(3)*exp(-k0(4)*(NPI(t,3)))) < 1
            r(t,ss) = 1/(1+exp(-k0(5)*(t-k0(6))));
        else
            r(t,ss) = (1/(k0(3)*exp(-k0(4)*(NPI(t,3)))))/(1+exp(-k0(5)*(t-k0(6))));
        end
    end
    I(1) = k0(8);
    E(1) = k0(7);
    S(1) = pop-E(1)-I(1)-vac(1);
    for t = 2 : tspan
        S(t) = S(t-1) - b(t-1,ss)*I(t-1)*(1-Controls(t-1))*S(t-1)/pop;
        E(t) = E(t-1) + b(t-1,ss)*I(t-1)*(1-Controls(t-1))*S(t-1)/pop - l*E(t-1);
        I(t) = I(t-1) + l*E(t-1) - r(t-1,ss)*I(t-1);
        I_det(t,ss) = r(t-1,ss)*I(t-1);
        I_new(t,ss) = l*E(t-1);
    end

end

%%
varNames  = {'Mean','CI05','CI95','duration','D_CI05','D_CI95','Sum','Sum_CI05','Sum_CI95','max','max_CI05','max_CI95',...
    'citycode','original_start','R0','Start.NPI','BI','FI','SI','C','Latent'};
varTypes = {'double','double','double','double','double','double','double','double','double','double','double','double',...
    'string','string','string','string','string','string','string','string','string'};
rowcol = [1,21];
sim = table('Size',rowcol ,'VariableTypes',varTypes,'VariableNames',varNames);

Inew = mean(I_new(2:end,:),2);
I05 = prctile(I_new(2:end,:),0.25,2);
I95 = prctile(I_new(2:end,:),99.75,2);
if isempty(I05(I05>=1))
    sim(1,2) = array2table(0);
    sim(1,5) = array2table(0);
    sim(1,8) = array2table(0);
    sim(1,11) = array2table(0);
else
    sim(1,2) = array2table(mean(I05(I05>=1)));
    sim(1,5) = array2table(max(I05(I05>=1)));
    sim(1,8) = array2table(sum(I05(I05>=1)));
    sim(1,11) = array2table(max(I05(I05>=1)));
end

sim(1,1) = array2table(mean(Inew(Inew>=1)));
sim(1,3) = array2table(mean(I95(I95>=1)));
sim(1,4) = array2table(max(Inew(Inew>=1)));
sim(1,6) = array2table(max(I95(I95>=1)));
sim(1,7) = array2table(sum(Inew(Inew>=1)));
sim(1,9) = array2table(sum(I95(I95>=1)));
sim(1,10) = array2table(max(Inew(Inew>=1)));
sim(1,12) = array2table(max(I95(I95>=1)));
sim(1,13) = {data.citycode(1)};
sim(1,14) = {convertTo(data.original_start(1),'yyyymmdd')};

%writetable(sim,[outdir,char('/Scenario_NPI=0'+string(data.citycode(1))+'.xlsx')],'WriteRowNames',true)




