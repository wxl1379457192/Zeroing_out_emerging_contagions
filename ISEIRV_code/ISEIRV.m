function [Effect, total,precases,valcases, RMSE_val,NRMSE_val, R2_val] = ISEIRV(data, NPIs, number, outdir)

% parameter initialization
pop = data.pop(1)*1e4;        % population

NPIs = smoothdata([NPIs(1,:); NPIs],'gaussian',5);

cases = data.I;
Control = data.Control;
Control = smoothdata(Control,'gaussian',14);

vac = pop*data.("Practical.vaccination");
vac(2:end) = vac(2:end) - vac(1:end-1);
%%
switch string(data.VG(1))
    case 'original&alpha'
        R0 = 4.475;
        l = 1/3; % 1/ latent period
        r0 = 7;
    case 'delta'
        R0 = 4.9;
        l = 1/2.7; % 1/ latent period
        r0 = 6;
    case 'omicron'
        R0 = 9.5;
        l = 1/2.3; % 1/ latent period
        r0 = 5;
end
b0 = R0/r0;
% sampling
b0 = makedist('normal','mu',b0,'sigma',1);
b1 = makedist('uniform','lower',0,'upper',3);%3
r0 = makedist('normal','mu',r0,'sigma',1);%1
r1 = makedist('uniform','lower',0,'upper',3);%3
r2 = makedist('uniform','lower',0,'upper',1);%1
r3 = makedist('uniform','lower',0,'upper',length(cases));
E0 = makedist('uniform','lower',3,'upper',15);%3-15
I0 = makedist('uniform','lower',1,'upper',5);%1-15

rate = 1;
variable_num = 8;
% training
train_step = 0;
train_num = 15;%15
steps = 10000;
samples = steps/5;
while (train_step<=train_num)
    MSE = zeros(steps,variable_num+1);
    for i = 1 : steps
        k0 = [b0.random b1.random r0.random r1.random r2.random r3.random E0.random I0.random];
        tspan = floor(length(cases)*rate)+1;
        if length(cases) - tspan < 4
           tspan = length(cases)-3;%3
        end
        r = zeros(tspan,1);
        b = zeros(tspan,1);
        for t = 1 : tspan
            b(t) = k0(1)*exp(-k0(2)*NPIs(t,1))*(1-0.25*NPIs(t,2));
            if (k0(3)*exp(-k0(4)*(NPIs(t,3)))) < 1
                r(t) = 1/(1+exp(-k0(5)*(t-k0(6))));
            else
                r(t) = (1/(k0(3)*exp(-k0(4)*(NPIs(t,3)))))/(1+exp(-k0(5)*(t-k0(6))));
            end
        end
        S = zeros(tspan,1);
        E = zeros(tspan,1);
        I = zeros(tspan,1);
        I_new = zeros(tspan,1);
        I_det = zeros(tspan,1);
        I(1) = k0(8);
        E(1) = k0(7);
        S(1) = pop-E(1)-I(1)-vac(1);
        for t = 2 : tspan
            S(t) = S(t-1) - b(t-1)*I(t-1)*(1-Control(t-1))*S(t-1)/pop;
            %S(t) = S(t-1) - b(t-1)*I(t-1)*S(t-1)/pop;
            E(t) = E(t-1) + b(t-1)*I(t-1)*(1-Control(t-1))*S(t-1)/pop - l*E(t-1);
            %E(t) = E(t-1) + b(t-1)*I(t-1)*S(t-1)/pop - l*E(t-1);
            I(t) = I(t-1) + l*E(t-1) - r(t-1)*I(t-1);
            I_new(t) = l*E(t-1);
            I_det(t) = r(t-1)*I(t-1);
        end
        MSE(i,:) = [k0 immse(I_new(2:tspan),cases(2:tspan))];
    end
    [~, I] = sort(MSE(:,end));
    
    % update the pdf
    k = MSE(I(1:samples),1:end-1);
    
    b0 = fitdist(k(:,1),'normal');
    b1 = fitdist(k(:,2),'normal');
    r0 = fitdist(k(:,3),'normal');
    r1 = fitdist(k(:,4),'normal');
    r2 = fitdist(k(:,5),'normal');
    r3 = fitdist(k(:,6),'normal');
    E0 = fitdist(k(:,7),'normal');
    I0 = fitdist(k(:,8),'normal');
    
    train_step = train_step + 1;
end


%% output
samples = 4000;
tspan = length(cases)+1;
I_new = zeros(tspan,samples);
I_det = zeros(tspan,samples);
S = zeros(tspan,1);
E = zeros(tspan,1);
I = zeros(tspan,1);
r = zeros(tspan,samples);
b = zeros(tspan,samples);
%MSE_val = zeros(1,samples);
%R2_val = zeros(1,samples);
Effect = zeros(samples,variable_num);
for ss = 1 : samples
    k0 = [b0.random b1.random r0.random r1.random r2.random r3.random E0.random I0.random];
    %
    Effect(ss,:) = k0;
    %
    tspan = length(cases)+1;
    for t = 1 : tspan
        b(t,ss) = k0(1)*exp(-k0(2)*NPIs(t,1))*(1-0.25*NPIs(t,2));
        if (k0(3)*exp(-k0(4)*(NPIs(t,3)))) < 1
            r(t,ss) = 1/(1+exp(-k0(5)*(t-k0(6))));
        else
            r(t,ss) = (1/(k0(3)*exp(-k0(4)*(NPIs(t,3)))))/(1+exp(-k0(5)*(t-k0(6))));
        end
    end
    I(1) = k0(8);
    E(1) = k0(7);
    S(1) = pop-E(1)-I(1)-vac(1);
    for t = 2 : tspan
        S(t) = S(t-1) - b(t-1,ss)*I(t-1)*(1-Control(t-1))*S(t-1)/pop;
        E(t) = E(t-1) + b(t-1,ss)*I(t-1)*(1-Control(t-1))*S(t-1)/pop - l*E(t-1);
        %S(t) = S(t-1) - b(t-1,ss)*I(t-1)*S(t-1)/pop;
        %E(t) = E(t-1) + b(t-1,ss)*I(t-1)*S(t-1)/pop - l*E(t-1);
        I(t) = I(t-1) + l*E(t-1) - r(t-1,ss)*I(t-1);
        I_det(t,ss) = r(t-1,ss)*I(t-1);
        I_new(t,ss) = l*E(t-1);
    end
    %Median_Inew = prctile(I_new(round(length(cases)*rate)+2:end,ss),50,2);
    %MSE_val(1,ss) = immse(Median_Inew,cases((round(length(cases)*rate)+1):end));%/sum(cases)
    %f = corrcoef(Median_Inew,cases((round(length(cases)*rate)+1):end));
    %R2_val(1,ss) = f(1,2,1);
end

if length(cases) - round(length(cases)*rate)+1 <4
    Median_Inew = prctile(I_new(end-4:end,:),50,2);
    valcases = cases(end-4:end);
    RMSE_val = sqrt(immse(Median_Inew,valcases));
    NRMSE_val = sqrt(immse(Median_Inew,valcases))/(sum(cases)/length(cases));
    f = corrcoef(Median_Inew,valcases);
    R2_val = f(1,2,1);
else
    Median_Inew = prctile(I_new(round(length(cases)*rate)+2:end,:),50,2);
    valcases = cases((round(length(cases)*rate)+1):end);
    RMSE_val = sqrt(immse(Median_Inew,valcases));
    NRMSE_val = sqrt(immse(Median_Inew,valcases))/(sum(cases)/length(cases));
    f = corrcoef(Median_Inew,valcases);
    R2_val = f(1,2,1);
end

precases = mean(sum(I_new));
valcases = sum(cases);
total = [mean(sum(I_new)), sum(cases)];
%% plot
tspan = length(cases)+1;
subplot(4,2,[1 2 3 4]);
p=imagesc(NPIs',[0,1]);
p.AlphaData=0.45;
yticks([1 2 3])
yticklabels({'Contact' 'Mask' 'Detection'})
colormap(summer)
colorbar('SouthOutside');
str = string(data.citycode(1))+'  From   '+string(data.Date(1))+'  to   '+string(data.Date(end))+' '+string(data.VG(1));
title(str,'FontName','times','FontSize',24)
yyaxis right
I_n_CI = prctile(I_new(2:end,:),[0.25 50 99.75],2);
I_n_Conf = [I_n_CI(:,1); I_n_CI(end:-1:1,3)];
p1=fill([1:tspan-1 tspan-1:-1:1],I_n_Conf','r');
p1.FaceColor = [1 0.8 0.8];
p1.EdgeColor = 'none';
p1.FaceAlpha = 0.5;
hold on
p3 = plot(1:tspan-1,cases,'b','LineWidth',1,'LineStyle','-.');
p2 = plot(1:tspan-1,mean(I_new(2:end,:),2),'b','LineWidth',1.5,'LineStyle','-');
legend([p2 p3],[{char('Estimation '+string(total(1)))}...
    {char('Infection '+string(total(2)))}],'Box','off','FontSize',16)

subplot(4,2,[5 6])
plot(mean(b,2),'LineWidth',1.5)
yyaxis right
plot(mean(r,2),'LineWidth',1.5)
legend([{'Transmission rate'} {'Detection rate'}],'location','best','FontSize',12)

subplot(4,2,7)
xspan = linspace(1,50,1000);
plot(xspan,E0.pdf(xspan),'LineWidth',1.5)
hold on
plot(xspan,I0.pdf(xspan),'LineWidth',1.5)
legend([{'Exposed at day 0'} {'Infections at day 0'}])

subplot(4,2,8)
xspan = linspace(0,10,1000);
plot(xspan,b0.pdf(xspan),'LineWidth',1.5)
hold on
plot(xspan,r0.pdf(xspan),'LineWidth',1.5)
legend([{'Basic transmission rate'} {'Basic infectious period'}])

set(gcf,'Position',[100 100 1500 1500])
%%
str = char('/Figs'+string(data.VG(1))+string(number));
saveas(gcf,[outdir, str],'png')
close()
vallength = length(cases) - round(length(cases)*rate);
cases = array2table(cases);
cases(:,2) = array2table(mean(I_new(2:end,:),2));
cases(:,3) = {string(data.VG(1))};
cases(:,4) = {data.citycode(1)};
cases(:,5) = {convertTo(data.original_start(1),'yyyymmdd')};
cases(:,6) = {vallength};
cases.Properties.VariableNames(1:6) = {'Observed_cases','Predicted_cases','VG','citycode','original_start','validation_length'};
writetable(cases,[outdir,char('/SEIR_simulation/'+string(data.VG(1))+string(data.citycode(1))+'.xlsx')],'WriteRowNames',true)

