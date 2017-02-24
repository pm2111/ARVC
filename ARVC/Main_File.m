close all; clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main File to run the O'Hara-Rudy Human Ventricular Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main File to run the O'Hara-Rudy human ventricular model. 
% Simulations results are saved as .mat in "Sim_Data", while figures are 
% saved in the "Sim_Figures\sim_name" folder, defined by the user.
%
% Simulations are run beat by beat. Only the last beat is saved
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Settings: 
% simulation name: change it for every different simulation
sim_name = 'OHR_endo_voltage_clamp_Na_condunctancex2';
% model_list: all the models that need to be tested/compared
% model_list = {'mod_1','mod_2','mod_3'...};
% where mod_n is the name of the function file containing the model
model_list = { 'modORd_ENDO_ARVC'};
% How many models? -> n_mod
n_mod = size(model_list,2);
% overwrite: 1 overwrite, 0 otherwise 
% (one for each model or a single value for all the models)
v_ow = [1];
% nb -> beats to run for each simulation:
% (one for each model or a single value for all the models)
v_nb = [150];            
% BCL -> basic cycle length for each simulation and modlen(ap_range)el (ms): 
% (one for each model or a single value for all the models)
v_BCL = [1000];  
%drug concentrations for testing, need stable APs up to 30x dose (0.1 microM)
drug_concentrations = 1; %[0.001 0.01 0.03 0.1 0.3 ]; %in micro M 
drug_name = 'amiodarone';
h_coef = [1.09 0.69 0.97];  % [Kr CaL Na ] put in this order 
IC50 = [0.86 1.9 15.9]; %[Kr CaL Na ] put in this order!!!!

n_con = size(drug_concentrations,2);
% data_save: 1 save Data and Figures, 0 otherwise
data_save = 0;
% Figure Flag (1 to produce the figure, 0 otherwise)
Fig_AP = 0;
Fig_ci =0;
Fig_I = 0;
fig_apd_response =0;
fig_drugs =0;
fig_ap_shape =0;
% Extracellular Concentrations
cNao = 140;
cCao = 1.8;
cKo = 5.4;
% ODE settings
ODEstep = 1; % ms 
options=odeset('MaxStep',ODEstep);
% Others settings
addpath(genpath('Models'));
addpath(genpath('Functions'));
legend_label = cell(n_mod,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Conditions
    v=-87;      nai=7;      nass=nai;       ki=150;         kss=ki;
    cai=1.0e-4; cass=cai;   cansr=1.2;      cajsr=cansr;    m=0;
    hf=1;       hs=1;       j=1;            hsp=1;          jp=1;
    mL=0;       hL=1;       hLp=1;          a=0;            iF=1;
    iS=1;       ap=0;       iFp=1;          iSp=1;          d=0;
    ff=1;       fs=1;       fcaf=1;         fcas=1;         jca=1;
    nca=0;      ffp=1;      fcafp=1;        xrf=0;          xrs=0;
    xs1=0;      xs2=0;      xk1=1;          Jrelnp=0;       Jrelp=0;
    CaMKt=0;
    
    y0_0 = [v       nai     nass    ki      kss...
            cai     cass    cansr   cajsr   m... 
            hf      hs      j       hsp     jp... 
            mL      hL      hLp     a       iF...
            iS      ap      iFp     iSp     d...
            ff      fs      fcaf    fcas    jca...
            nca     ffp     fcafp   xrf     xrs...
            xs1     xs2     xk1     Jrelnp  Jrelp...
            CaMKt]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Directories to save/load simulation results and figures:
% Simulation Data
dir_Sim = 'Sim_Data';
% Figures
dir_Fig = 'Sim_Figures';
% Save
dir_Save = [dir_Fig '/' sim_name];
% If it exists, ask to overwrite
if ~exist(dir_Save,'dir')
 mkdir(dir_Save); flag=1;
else
 flag=input('The Sim Directory exists, overwrite?\n 1 -> Yes, 0 -> No:\n');
end
% FLAG
if flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise cells that should not be overwritten in simu step 
APDs = cell(n_mod,n_con); %can also include different apd fractions (not included here)
APD_set = [30 50 70 90];
%save all AP curves for plotting 
APs = cell(n_mod);
ts = cell(n_mod);
voltage_response=cell(2);

%% Simulations
for v_cl = 1:30 %iterate models over a range of constant voltage stimuli
             voltage_response{1,v_cl} = -100+5*v_cl; %store voltage clamps used 
             y0_0(1) = voltage_response{1,v_cl};
for i_mod = 1:n_mod
    for j_mod = 1 :n_con
        disp(['Sim model number ',num2str(i_mod),'out of ', n_mod,model_list{i_mod},'stim', num2str(y0_0(1)), ' mV']);
        mod =  str2func(model_list{i_mod});
        if length(v_ow)==1; ow=v_ow; else ow=v_ow(i_mod); end;
        if length(v_nb)==1; nb=v_nb; else nb=v_nb(i_mod); end;
        if length(v_BCL)==1; BCL=v_BCL; else BCL=v_BCL(i_mod); end;
        concen = drug_concentrations(j_mod); %select model drug concentration 
        % Data loading
        dir_Sim_Save = [dir_Sim '/' model_list{i_mod}];
        if ~exist(dir_Sim_Save,'dir')
            mkdir(dir_Sim_Save);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Figures Inizialization (first model only):
        if i_mod == 1 
           % Figure 1: Action Potential
           if Fig_AP == 1
            hAP = figure('Name','Action Potential'); hold on; xlabel('t (ms)'); ylabel('V_m (mV)'); box;
           end
           % Figure 2: Intracellular Concentrations
           if Fig_ci == 1 && nb>1
            hci = figure('Name','Intracellular Concentrations'); 
            subplot(3,1,1); hold on; xlabel('beats'); ylabel('[Na]_i'); box;
            subplot(3,1,2); hold on; xlabel('beats'); ylabel('[Ca]_i'); box;
            subplot(3,1,3); hold on; xlabel('beats'); ylabel('[K]_i'); box; 
           end
           if Fig_I==1
           % Figure 3: Ionic Current I (Istim, INa, INaL & ICaL)
            hI1 = figure('Name','Ionic Currents - I');
            subplot(2,2,1); hold on; xlabel('t (ms)'); ylabel('I_{Stim} pA/pF'); box; 
            subplot(2,2,2); hold on; xlabel('t (ms)'); ylabel('I_{Na} pA/pF'); box; 
            subplot(2,2,3); hold on; xlabel('t (ms)'); ylabel('I_{CaL} pA/pF'); box; 
            subplot(2,2,4); hold on; xlabel('t (ms)'); ylabel('I_{NaL} pA/pF'); box; 
           % Figure 4: Ionic Current II (Ito, IK1, IKr, IKs)
            hI2 = figure('Name','Ionic Currents - II');
            subplot(2,2,1); hold on; xlabel('t (ms)'); ylabel('I_{to} pA/pF'); box; 
            subplot(2,2,2); hold on; xlabel('t (ms)'); ylabel('I_{Kr} pA/pF'); box; 
            subplot(2,2,3); hold on; xlabel('t (ms)'); ylabel('I_{K1} pA/pF'); box; 
            subplot(2,2,4); hold on; xlabel('t (ms)'); ylabel('I_{Ks} pA/pF'); box; 
           % Figure 5: Ionic Current III (INaK, INaCa, Jrel, Jup)
            hI3 = figure('Name','Ionic Currents - III');
            subplot(2,2,1); hold on; xlabel('t (ms)'); ylabel('I_{NaK} pA/pF'); box; 
            subplot(2,2,2); hold on; xlabel('t (ms)'); ylabel('J_{rel} mM/ms'); box; 
            subplot(2,2,3); hold on; xlabel('t (ms)'); ylabel('I_{NaCa} pA/pF'); box; 
            subplot(2,2,4); hold on; xlabel('t (ms)'); ylabel('J_{up} mM/ms'); box; 
            % Figure 6: Ionic Current III (INaK, INaCa, Jrel, Jup)
            hI4=  figure('Name', 'Short/Long Channels');        
            subplot(2,2,1); hold on; xlabel('t (ms)'); ylabel('I_{NaCa i} pA/pF'); box; 
            subplot(2,2,2); hold on; xlabel('t (ms)'); ylabel('I_{Na peak} pA/pF'); box; 
            subplot(2,2,3); hold on; xlabel('t (ms)'); ylabel('I_{Ca peak} pA/pF'); box; 
            %subplot(2,2,4); hold on; xlabel('t (ms)'); ylabel('I_{NaL} pA/pF'); box; 
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vNai=zeros(1,nb);  
        vCai=zeros(1,nb);   
        vKi=zeros(1,nb);   
        % model inputs
        input_args ={};
        input_args{1} = 1;
        input_args{2} = 2;
        input_args{3} = BCL;
        input_args{4} = [cNao cCao cKo]; 
        input_args{5} = concen;
        input_args{6} = h_coef;
        input_args{7} =IC50;
        % File to save Sim Data
        data_file = [model_list{i_mod} '_' num2str(nb) 'b_' num2str(BCL) drug_name 'concen' num2str(drug_concentrations(j_mod)) 'ms.mat'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Simulations
        % If the simulation data exist, load them:
        if exist([dir_Sim_Save '/' data_file],'file') && ow==0
            load([dir_Sim_Save '/' data_file],'t','y','CVs','intra_c');
            fprintf('Simulation Data loaded\n');
            vNai=intra_c(1,:);
            vCai=intra_c(2,:);
            vKi=intra_c(3,:);
        % If not, run a new simulation:
        else
            y0 = y0_0;
        % all beats are simulated for BCL
        
        for n = 1:nb
            [t,y] = ode15s(mod,[0 BCL],y0,options,input_args{:});
            %fprintf('Beat %i of %i\n',n,nb);
            if n ==nb
            APs{i_mod,j_mod} = y(:,1); %output last beat AP for every drug concen, last one has no drug
            ts{i_mod, j_mod} = t;
            end
            %APDsname =  ['APD Duration at' num2str(drug_concentrations(i)) 'microM'];

            APD = computeAPD_rudy(t,y(:,1),90,0);
            APDs{i_mod,j_mod,n} = APD;
        
            % new initial conditions:
            y0 = y(end,:);   
            vNai(n)=y(end,2);  
            vCai(n)=y(end,6); 
            vKi(n)=y(end,4);  
           
        end


        % Computed Variables (e.g. ionic currents)
        input_args{1}=0; 
        lCVs=size(feval(mod,t(1),y(1,:),input_args{:}),1);
        CVs = zeros(length(t),lCVs);
        for j=1:size(y,1)
            CVs(j,:)=feval(mod,t(j),y(j,:),input_args{:});
        end
        intra_c = [vNai;vCai;vKi];
        save([dir_Sim_Save '/' data_file],'t','y','CVs','intra_c');    
        end
        
        
        
        voltage_response{2,v_cl} = CVs(end,2); %only sodium current is studied
        disp(CVs(end,2));
        end %ends voltage clamp loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Figures Plot:
        if Fig_AP == 1
            figure(hAP); 
            plot(t,y(:,1),'LineWidth',2); 
        end
        % Figure 2: Intracellular Concentrations
        if Fig_ci == 1 && nb>1
            figure(hci);  
            subplot(3,1,1); plot(1:nb,vNai,'LineWidth',2);
            subplot(3,1,2); plot(1:nb,vCai,'LineWidth',2);
            subplot(3,1,3); plot(1:nb,vKi,'LineWidth',2);
        end
        if Fig_I==1
            % Figure 3: Ionic Current I (Istim, INa, INaL & ICaL)
            figure(hI1)
            subplot(2,2,1); plot(t,CVs(:,1),'LineWidth',2);
            subplot(2,2,2); plot(t,CVs(:,2),'LineWidth',2);
            subplot(2,2,3); plot(t,CVs(:,8),'LineWidth',2);
            subplot(2,2,4); plot(t,CVs(:,3),'LineWidth',2);
            % Figure 4: Ionic Current II (Ito, IK1, IKr, IKs)
            figure(hI2)
            subplot(2,2,1); plot(t,CVs(:,4),'LineWidth',2);
            subplot(2,2,2); plot(t,CVs(:,5),'LineWidth',2);
            subplot(2,2,3); plot(t,CVs(:,7),'LineWidth',2);
            subplot(2,2,4); plot(t,CVs(:,6),'LineWidth',2);
            % Figure 5: Ionic Current III (INaK, INaCa, Jrel, Jup)
            figure(hI3)
            subplot(2,2,1); plot(t,CVs(:,9),'LineWidth',2);
            subplot(2,2,2); plot(t,CVs(:,11),'LineWidth',2);
            subplot(2,2,3); plot(t,CVs(:,10),'LineWidth',2);
            subplot(2,2,4); plot(t,CVs(:,12),'LineWidth',2);
            %Figure 6: Late and peak Na and Ca
            figure(hI4)
            subplot(2,2,1); plot(t,CVs(:,13),'LineWidth',2);
            subplot(2,2,2); plot(t,CVs(:,16),'LineWidth',2);
            subplot(2,2,3); plot(t,CVs(:,15),'LineWidth',2);

        end
        legend_label{i_mod} = [char(model_list{i_mod}) ' ' num2str(nb) 'b ' num2str(BCL) 'ms'];


        ap_range = [1 15 150 500 600 700 800 900 1000];
        cols=hsv(max(size(ap_range))); %create a diverse colour palette
    %     figure('Name','AP stability')
    %     hold all;
    %     for z = 1:max(size(ap_range))
    %          plot(t(1:1000),APs{1,ap_range(z)}(1:1000), 'Linewidth',2,'Color',cols(z,:) , 'DisplayName', strcat('beat number', num2str(ap_range(z))));
    %     end
    %     grid on;function output=modORd_ENDO(t,X,varargin) 

%compute APD averages

for i = 1:i_mod
    for j = 1:4
        array = [APDs{i,j_mod,j,150:end}];
        av_APD(i,j) = mean(array);
    end 
end
%curve fit 
%fitobject = fit(transpose(a(1:4)), transpose(av_APD(1,:)),'exp1');
%APD figure in O Hara Rudy paper (fig #7)
apd_color = ['y', 'm','c', 'r', 'g','b','w', 'k'];
if fig_apd_response ==1
figure('Name','APD response to CL')
hold all;
for i = 1:i_mod
    for j =1:max(size(APD_set))
        
    scatter(v_BCL(i), av_APD(i,j),apd_color(j), 'Linewidth',2,'DisplayName', strcat('APD %', num2str(APD_set(j))));
    end
end
legend('show');
end
if fig_drugs == 1
    for j = 1:n_con
        name =  ['APD Duration for ' drug_name ' at ' num2str(drug_concentrations(j)) 'microM'];
        figure('Name', name)
        hold all;
        
        for i=1:n_mod
                if i < n_mod
                    plot([APDs{i,j,:}],apd_color(i), 'Linewidth',2,'DisplayName', strcat('CL', num2str(v_BCL(i))));
                
                else
                    plot([APDs{i,j,:}],'bl', 'Linewidth',4,'DisplayName', strcat('no drug, CL', num2str(1000)));
                end
        end
         
        xlabel('beat number');
        ylabel('APD90 Duration (ms)');
        grid on;
        legend('show');
        %saveas(name,[dir_Save '/' name]);
    end
end
  if fig_ap_shape ==1
    for j = 1:n_con
    figure('Name', ['AP shape for ' drug_name ' at ' num2str(drug_concentrations(j)) 'microM'])
    hold all;
    for i=1:n_mod
        if i < n_mod
                plot(ts{i,j,:},[APs{i,j,:}],apd_color(i), 'Linewidth',2,'DisplayName', strcat('CL', num2str(v_BCL(i))));
        else
                plot(ts{i,j,:},[APs{i,j,:}],'b','--', 'Linewidth',4,'DisplayName', strcat('no drug, CL', num2str(1000)));
        end
    end
        xlabel('time(ms)');
        ylabel('Membrane potential (mV)');
        grid on;
        legend('show');
        %saveas(name,[dir_Save '/' name]);
    end
        
 end
        
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure Legend & Saving
if Fig_AP==1
    figure(hAP); legend(char(legend_label)); saveas(hAP,[dir_Save '/APtraces']);
end
if Fig_ci==1 && nb>1
    figure(hci); legend(char(legend_label)); saveas(hci,[dir_Save '/IntraConc']);
end
if Fig_I==1
    figure(hI1); legend(char(legend_label)); saveas(hI1,[dir_Save '/IonicCurrents_1']);
    figure(hI2); legend(char(legend_label)); saveas(hI1,[dir_Save '/IonicCurrents_2']);
    figure(hI3); legend(char(legend_label)); saveas(hI1,[dir_Save '/IonicCurrents_3']);
end   

%%%%%%%%%%%%%%%
%%%SAVING DATA POINTS
dlmwrite(sim_name, voltage_response);
v100 = csvread('Sim_Data/endo_Na_V_clamp');
figure('Name', 'Na Conductance Dependent Current')
hold all;
plot([voltage_response{1,1:end}], [voltage_response{2,1:end}],'LineWidth',2.0, 'DisplayName', 'Na Conductance reduced by 36%');
plot(v100(1,:), v100(2,:),'LineWidth',2.0, 'DisplayName', 'Healthy Na Conductance');
legend('show');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% END FLAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%