function [diamString] = calcEvap(diam_0,u_slip,t_evap,strSize)
%Adriatique - Craft
%Initialization and control parameters
prop_opt  = 2;
model_VLE = 3;
diamString = zeros(strSize,1);

gas_input           = 1;        %Vary between air and nitrogen

input_tb            = prop_opt;       
input_liqdensity    = prop_opt;
input_liqlambda     = prop_opt;
input_bindiff       = prop_opt; 
input_gascp         = prop_opt; %Be aware when changing gas_input
input_vaporcp       = prop_opt;
input_vaporlambda   = prop_opt; 
input_gaslambda     = prop_opt;
input_visc          = prop_opt;
input_lv            = prop_opt;
input_cl            = prop_opt;
nu_sh_opt           = prop_opt;

mat_decane          = 1;
mat_ethanol         = 2;
mat_water           = 3;
mat_hexane          = 4;
mat_heptane         = 5;

gas_air             = 1;
gas_N2              = 2;


PG       = 101325.0e0;
Patm     = 101325.0e0;

% Gas settings
if gas_input==1
    WG       = 28.97e0;
elseif gas_input==2
    WG       = 28.0134e0;
end    

TG       = 294.0e0;
Td       = 307.0e0;
YG       = 0.008541;
diamp    = diam_0*1e-6;
subst    = mat_water;

WV       = Mw(subst);
theta2   = WG/WV;
R        = 8314.3e0;

Tb       = T_boiling(subst,input_tb);
T_wb     = 137*(Tb/373.15)^0.68*log10(TG)-45;
rule     = 1.0d0/3.0d0;
heat_up  = Td/T_wb; %if heat_up~=1, then we have a straight line (D^2 law)

dt       = 0.001;
tfinal   = t_evap(end);
i        = 1;
count_t  = 1;
ttotal   = 0.0;
evap_tol = 1.0d-10;
tol      = 10*evap_tol;
maxTvar  = 10d0;

%Step 0 - Initializing rho_l and md
rho_l            = liqDensity(subst,Td,input_liqdensity);
md               = pi/6.0d0*rho_l*diamp^3.0d0;

%These steps are performed here before the loop to compute the initial u_slip
%Step 1 - Vapor Pressure and Mass/Molar fractions
LV               = heat_vaporization(subst,Td,input_lv);
vappres          = vapor_pressure(subst,Td,model_VLE,input_tb,LV);
Xseq             = vappres/max(PG,1.0d-15);
Xseq             = min(max(Xseq,1d-15),1.0d0-1d-15);
Yseq             = Xseq/(Xseq+(1.0d0-Xseq)*theta2);
Yseq             = min(max(Yseq,1d-15),1.0d0-1d-15);

%Step 2 - Film equivalent Mass/Molar Fractions and Temperature
YR               = film_rule(YG,Yseq,rule);
TR               = film_rule(TG,Td,rule);
XR               = (theta2)*YR/(1.0e0+(theta2-1.0e0)*YR);

%Intermediate Step - Computing u_slip (supposed constant)
rho_f            = film_density(TR,YR,PG,WV,WG);
[mu_f,~,~]       = film_visc(subst,XR,WV,WG,TR,input_visc,gas_input);

Re_d             = u_slip*rho_f*diamp/mu_f;
%% 

while(ttotal<tfinal&&diamp>2e-6)
    
    %Step 1 - Vapor Pressure and Mass/Molar fractions
    LV               = heat_vaporization(subst,Td,input_lv);
    vappres          = vapor_pressure(subst,Td,model_VLE,input_tb,LV);
    Xseq             = vappres/max(PG,1.0d-15);
    Xseq             = min(max(Xseq,1d-15),1.0d0-1d-15);
    Yseq             = Xseq/(Xseq+(1.0d0-Xseq)*theta2);
    Yseq             = min(max(Yseq,1d-15),1.0d0-1d-15);
    
    %Step 2 - Film equivalent Mass/Molar Fractions and Temperature
    YR               = film_rule(YG,Yseq,rule);
    XR               = (theta2)*YR/(1.0e0+(theta2-1.0e0)*YR);
    TR               = film_rule(TG,Td,rule);
    
    %Step 3 - Film properties
    Gamma_f          = film_bindiff(subst,TR,WV,WG,Patm,input_bindiff,gas_input);
    [mu_f,~,~]       = film_visc(subst,XR,WV,WG,TR,input_visc,gas_input);
    rho_f            = film_density(TR,YR,PG,WV,WG);
    [lambda_f,~,~]   = film_lambda(subst,XR,WV,WG,TR,input_vaporlambda,gas_input,input_gaslambda);
    [~,~,cp_f]       = film_cp(subst,TR,YR,input_vaporcp,gas_input,input_gascp);
    Sc_f             = mu_f/(Gamma_f*rho_f);
    Pr_f             = mu_f*cp_f/lambda_f;
    Le_f             = Sc_f/Pr_f;
    
    %Step 4 - Gas properties
    [~,~,mu_G]       = film_visc(subst,XR,WV,WG,TR,input_visc,gas_input);
    rho_G            = PG*WG/((R)*TR+1.0e-15);
    [~,~,lambda_G]   = film_lambda(subst,XR,WV,WG,TR,input_vaporlambda,gas_input,input_gaslambda);
    [~,cp_G,~]       = film_cp(subst,TR,YR,input_vaporcp,gas_input,input_gascp);
    Sc_G             = mu_G/(Gamma_f*rho_G);    
    Pr_G             = mu_G*cp_G/lambda_G;
    %Le_G             = Sc_G/Pr_G;
    
    %Step 5 - Liquid properties
    rho_l            = liqDensity(subst,Td,input_liqdensity);
    cl               = liquid_cp(subst,Td,input_cl);
    
    %Step 6 - Intermediate calculations
    Re_d             = rho_f*u_slip*diamp/mu_f; %mudou de G pra f
    cp_v             = vapor_cp(subst,TR,input_vaporcp);
    taup             = rho_l*diamp^2.0d0/(18.0d0*mu_f);
    theta1           = cp_v/(cl+1d-15); %testar cp_G ao inv�s de cp_v
    
    [Nu,Sh]          = nu_sh(Re_d,Pr_f,Sc_f,nu_sh_opt); %testar tamb�m Pr_G e Sc_G
    
    cp_bar           = cp_f; %Notice that cp_f=YR*cp_v+(1-YR)*cp_G already
    Le_bar           = Le_f;
    
    %Step 7 - Spalding Numbers and M2 Coefficients (H_M, H_DeltaT)
    BMeq             = (Yseq-YG)/(1.0d0-Yseq);
    BT               = (TG-Td)*cp_v/(max(LV,1.0d-15));
    H_M              = log(1.0d0+BMeq);
    
    %Step 8 - Spalding Number BT' loop
    BT_line          = BT;
    FM               = (1.0d0+BMeq)^0.7d0/BMeq*log(1.0d0+BMeq);
    Sh_star          = 2.0d0+(Sh-2.0d0)/(max(FM,1d-15));
    
    count=1;
    tol = 10.0;
    while(tol>1e-6&&count<10||count<=3)
        BT_old       = BT_line;
        FT           = (1.0d0+BT_line)^0.7d0/(BT_line+1d-15)*log(1.0d0+BT_line);
        Nu_star      = 2.0d0+(Nu-2.0d0)/(FT+1d-15);
        phi          = cp_v/cp_bar*Sh_star/Nu_star*1.0d0/Le_bar;
        BT_line      = (1.0d0+BMeq)^phi-1.0d0;
        tol          = abs((BT_line-BT_old)/(BT_old+1.0d-15));
        count        = count+1;
    end
    
    %Step 9 - Computing md_dot and updating md
    md_dot           = -Sh_star/max(3.0d0*Sc_f,1.0d-15)*(md/max(taup,1d-15))*H_M;
    md = md + md_dot*dt;
    md = max(md,1.0d-15);
    
    %Step 10 - Computing Td_dot and updating Td
    f2               = -md_dot/max(md*BT_line,1.0d-15)*(3.0d0*Pr_f*taup/max(Nu_star,1.0d-15));
    Td_dot           = f2*Nu_star/max(3.0d0*Pr_f,1d-15)*(theta1/max(taup,1d-15))*(TG-Td)+LV/max(cl,1d-15)*md_dot/max(md,1d-30);
    if (abs(Td_dot*dt)<=maxTvar)
        Td  = Td + Td_dot*dt;
        Td  = min(max(Td,Liq_Tmin(subst)),Tb-1.0d-3);
    else
        dt  = dt/2;
    end
    
    %Step 11 - Updating the Diameter
    rho_l   = liqDensity(subst,Td,input_liqdensity);
    diamp   = (6.0d0/(pi*rho_l)*md)^(1.0d0/3.0d0);
    
    Dp(i)   = diamp*1000;
    Tp(i)   = Td;
    
    t(i)    = ttotal;
    ttotal  = ttotal+dt;
    i       = i+1;
    
    dt_eRef = t_evap(2) - t_evap(1);
    
    if (ttotal <= t_evap(end) && count_t <= size(t_evap,2))
        if(ttotal>t_evap(count_t)-0.1*dt_eRef && ttotal<t_evap(count_t)+0.1*dt_eRef && marker == false)
%         if(round(ttotal,2)==round(t_evap(count_t),2))  
            diamString(count_t) = Dp(i-1);      
            count_t = count_t + 1;
            marker = true;
        else
            marker = false;
        end
    end
end

%Plots
% 
%     figure(1)
%     hold on
%     plot(t,Dp.*Dp)
%     title('Diameter squared (m^2) x Time (s)')
%     hold off
%     
%     figure(2)
%     hold on
%     plot(t,Tp)
%     title('Temperature (K) x Time (s)')
%     hold off
end
