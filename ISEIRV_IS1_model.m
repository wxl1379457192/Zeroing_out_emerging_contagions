function [Effect, total] = ISEIRV_sensitivity_analysis_IS1(data, NPIs)

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
            b(t) = k0(1)*exp(-k0(2)*NPIs(t,1))*(1-0.1*NPIs(t,2));
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
        b(t,ss) = k0(1)*exp(-k0(2)*NPIs(t,1))*(1-0.1*NPIs(t,2));
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

precases = mean(sum(I_new));
valcases = sum(cases);
total = [mean(sum(I_new)), sum(cases)];


