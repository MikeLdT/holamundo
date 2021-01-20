%% Code to reproduce the results in Appendix A of the paper 
%% "The Rise of Services and Balanced Growth in Theory and Data"
%% Miguel Leon-Ledesma and Alessio Moro
%% American Economic Journal Macroeconomics 2019

%% Initializing code

clear
clc

time = 66;

Time = 2015-time+1:2015;


%% Vectors for variables (see below for definitions)

e            = zeros(time,1);  
y            = zeros(time,1);  
k            = zeros(time,1);  
AI           = zeros(time,1);
As           = zeros(time,1);
ps           = zeros(time,1);
pI           = zeros(time,1);
pg          = zeros(time,1);
k            = zeros(time,1);
e            = zeros(time,1);
r            = zeros(time,1);
w            = zeros(time,1);
cg           = zeros(time,1);
cs           = zeros(time,1);
budget       = zeros(time,1);
expenditure  = zeros(time,1);
nom_GDP      = zeros(time,1);

%% Parameter Values


beta     = 0.95;                     % discount rate
alpha    = 0.34;                     % capital share in cobb-douglas
delta    = 0.06;                     % depreciation

nu      = 0.6350;
epsilon = 0.1699;
gamma   = 0.5031;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AI(1)   = 1;                        % initial TFP in investment
As(1)   = 1;                        % initial TFP in services
Ag(1)   = 1;                        % initial TFP in consumption goods
ps(1)   = (AI(1)/As(1))^(1-alpha);  % relative price services/manufacturing
pI(1)   = 1;                        % price of investment (numeraire in each period)   
pg(1)   = (AI(1)/Ag(1))^(1-alpha);  % price of consumption goods

gamma_I = 0.0242;                    % growth of tfp in the investment sector
gamma_s = gamma_I-0.0110/(1-alpha);  % growth of tfp in services. this is calibrated from gamma_I and the 
                                     % relative price services/investment in the data 
gamma_g =gamma_s+0.0161/(1-alpha);   % growth of tfp in goods. this is calibrated from gamma_s and the 
                                     % relative price services/goods in the data 

%% Initial capital   

k(1)     = (alpha*beta/((1+gamma_I)^(1-alpha*epsilon)*(1+gamma_s)^((alpha-1)*epsilon)-beta*(1-delta)))^(1/(1-alpha))*AI(1);                      

%% Equilibrium at time 1

e(1)        = AI(1)^(1-alpha)*k(1)^alpha-(gamma_I+delta)*k(1);          % total expenditure in goods
r(1)        = 1/beta*(1+gamma_I)^(1-alpha*epsilon)*(1+gamma_s)^((alpha-1)*epsilon)-1+delta; % BGP real interest rate in investment units
w(1)        = (1-alpha)*AI(1)^(1-alpha)*k(1)^alpha;                     % wage rate in invetment units
cg(1)       = nu*e(1)^(1-epsilon)*ps(1)^(epsilon-gamma);                % consumption of goods
cs(1)       = e(1)/ps(1)-nu*e(1)^(1-epsilon)*ps(1)^(epsilon-gamma-1);   % consumption of services
nom_GDP(1)  = AI(1)^(1-alpha)*k(1)^alpha;                               % "nominal" GDP (in investment units)

exp_ratio(1) = e(1)/nom_GDP(1);    % share of consumption in nominal GDP
invs(1)      = nom_GDP(1)-e(1);    % investment


%% Simulation for "time" periods


for t = 1 : time-1;
      
    AI(t+1)     = AI(t)*(1+ gamma_I); 
    As(t+1)     = As(t)*(1+ gamma_s); 
    Ag(t+1)     = Ag(t)*(1+ gamma_g); 
        
    ps(t+1)    = (AI(t+1)/As(t+1))^(1-alpha);
    pI(t+1)    = 1;
    pg(t+1)    = (AI(t+1)/Ag(t+1))^(1-alpha);
    
    k(t+1)     = k(t)*(1+ gamma_I); 
    e(t+1)     = e(t)*(1+ gamma_I); 
    
    r(t+1)     = 1/beta * (1+gamma_I)^(1-alpha*epsilon) * (1+gamma_s)^((alpha-1)*epsilon) - 1 + delta;
    
    w(t+1)     = (1-alpha)*AI(t+1)^(1-alpha)*k(t+1)^alpha;
    
    cg(t+1)    = nu*e(t+1)^(1-epsilon)*ps(t+1)^(epsilon-gamma)*pg(t+1)^(gamma-1);
    cs(t+1)    = e(t+1)/ps(t+1)-nu*e(t+1)^(1-epsilon)*pg(t+1)^gamma*ps(t+1)^(epsilon-gamma-1);
         
    budget (t+1) = e(t)+k(t+1)-(1-delta)*k(t)- AI(t)^(1-alpha)*k(t)^alpha; % variable to check equilibrium is verified. if = 0 ok
    expenditure_check(t+1) = e(t+1)-pg(t+1)*cg(t+1)-ps(t+1)*cs(t+1);       % variable to check equilibrium is verified. if = 0 ok
        
    nom_GDP(t+1)  = AI(t+1)^(1-alpha)*k(t+1)^alpha;
    
    invs(t+1)     = nom_GDP(t+1)-e(t+1);
    
    exp_ratio(t+1)= e(t+1)/nom_GDP(t+1);
    
    invs_out(t+1) = invs(t+1)/nom_GDP(t+1); 
    
end

%% Computing additional variables

GDP         = CompGDP2(time,ps,cs,pI,invs,pg,cg); % calls the function CompGDP2 which computes GDP from the model's outcome
Consumption = CompCons1(time,ps,cs,pg,cg);        % calls the function CompCons1 which computes Consumption from the model's outcome

price_GDP  = nom_GDP./(GDP);    % computes GDP deflator
price_Cons = e./(Consumption);  % computes Consumption deflator

real_rate      = r./price_GDP/(r(1)./price_GDP(1))  ; % computes the real interest rate in GDP units
real_rate_cons = r./price_Cons/(r(1)./price_Cons(1)); % computes the real interest rate in consumption units

MPKGDP = r./price_GDP/(r(1)./price_GDP(1))  ;   % computes the MPK in terms of GDP units 
MPKC   = r./price_Cons/(r(1)./price_Cons(1));   % computes the MPK in terms of Cons units 

k_gdp_ratio   = k./(GDP);                     % capital/output ratio. This is increasing, and consistent with Whelan 2002
invs_out_real = invs'./(GDP);                 % Investment/output ratio measured as Real Inv / Real GDP (i.e. using goods prices for I and GDP deflator for DDP)
GDP_growth    = GDP(2:time)./GDP(1:time-1)-1; % computes GDP growth

share = ps.*cs./e;                % computes the share of services in expenditure
share_GDP = ps.*cs./nom_GDP;      % computes the share of services in GDP



%% Data to plot

[data,text]=xlsread('data_to_plot','Foglio1');

%% Figure 10 in the paper

figure
subplot(1,3,1), plot(Time,log(GDP), 'b--',Time,log(data(:,1)), 'r','LineWidth',1.5)
title('log(GDP)')
set(gca,'fontsize',16)
axis([1950 2015 0 1.5])
subplot(1,3,3), plot(Time,log(invs_out_real./invs_out_real(1))+log(data(1,6)*100), 'b--', Time, log(data(:,6)*100) , 'r' ,'LineWidth',1.5)
legend('Model', 'Data')
title('log(Investment/Output Ratio)')
set(gca,'fontsize',16)
axis([1950 2015 2.6 3.3])
subplot(1,3,2), plot(Time,share , 'b--', Time,data(:,3) , 'r','LineWidth',1.5)
title('Services share in Cons.')
set(gca,'fontsize',16)
axis([1950 2015 0.35 0.7])


%% Data Statistics in Table 6

Goods_GDPpc_growth_data      = 0.0242;
Initial_services_share_data  = 0.393;
Final_services_share_data    = 0.6850;

Average_growth_real_investment_output_data  = 0.0030;
ps_pg_growth_data                           = 0.0161;
ps_pI_growth_data                           = 0.0110;


%% Model's Statistics in Table 6

Goods_GDPpc_growth     = gamma_I;
Initial_services_share = share(1);
Final_services_share   = share(66);

Average_growth_real_investment_output = (invs_out_real(66)/invs_out_real(1))^(1/(65))-1;
ps_pg_growth                          = (ps(time)./pg(time))^(1/(65))-1;
ps_pI_growth                          = (ps(time)./pI(time))^(1/(65))-1;


disp('Table 6 reports statistics for the data (row 1) and model (row 2) in the following order:' )
disp('Goods GDPpc growth; Initial share of services; Final share of services;' )
disp('Real I/Y Growth; Growth of ps/pg; Growth of ps/pI' )


Table_6 = [Goods_GDPpc_growth_data Initial_services_share_data Final_services_share_data...
           Average_growth_real_investment_output_data ps_pg_growth_data ps_pI_growth_data;
           Goods_GDPpc_growth Initial_services_share Final_services_share...
           Average_growth_real_investment_output ps_pg_growth ps_pI_growth]
       
Model_average_GDP_growth      = mean(GDP_growth)
Model_Nominal_Investment_Rate = (invs_out(2))

%% Decline in the MPK
MPK_GDP_decline          = (1-MPKGDP(66))
MPK_Consumption__decline = (1-MPKC(66))


