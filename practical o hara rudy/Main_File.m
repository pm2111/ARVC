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
sim_name = 'tempSIM';
% model_list: all the models that need to be tested/compared
% model_list = {'mod_1','mod_2','mod_3'...};
% where mod_n is the name of the function file containing the model
model_list = {'modORd_ENDO'};
% How many models? -> n_mod
n_mod = size(model_list,2);
% overwrite: 1 overwrite, 0 otherwise 
% (one for each model or a single value for all the models)
v_ow = [0];
% nb -> beats to run for each simulation:
% (one for each model or a single value for all the models)
v_nb = [15];            
% BCL -> basic cycle length for each simulation and model (ms): 
% (one for each model or a single value for all the models)
v_BCL = [1000];
% data_save: 1 save Data and Figures, 0 otherwise
data_save = 0;
% Figure Flag (1 to produce the figure, 0 otherwise)
Fig_AP = 1;
Fig_ci = 1;
Fig_I = 1;
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
if ~exist(dir_Save,'dir'); 
 mkdir(dir_Save); flag=1;
else
 flag=input('The Sim Directory exists, overwrite?\n 1 -> Yes, 0 -> No:\n');
end
% FLAG
if flag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulations
for i_mod = 1:n_mod
    fprintf('Sim %i of %i: %s\n',i_mod,n_mod,model_list{i_mod});
    mod =  str2func(model_list{i_mod});
    if length(v_ow)==1; ow=v_ow; else ow=v_ow(i_mod); end;
    if length(v_nb)==1; nb=v_nb; else nb=v_nb(i_mod); end;
    if length(v_BCL)==1; BCL=v_BCL; else BCL=v_BCL(i_mod); end;
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
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vNai=zeros(1,nb);  
    vCai=zeros(1,nb);   
    vKi=zeros(1,nb);   
    % model inputs
    input_args ={};
    input_args{1} = 1;
    input_args{2} = 1;
    input_args{3} = BCL;
    input_args{4} = [cNao cCao cKo]; 
    % File to save Sim Data
    data_file = [model_list{i_mod} '_' num2str(nb) 'b_' num2str(BCL) 'ms.mat'];
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
        fprintf('Beat %i of %i\n',n,nb);
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
    end
    legend_label{i_mod} = [char(model_list{i_mod}) ' ' num2str(nb) 'b ' num2str(BCL) 'ms'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% APD
    APD = computeAPD(t,y(:,1),90,0);
    fprintf('APD_90: %4.2f ms\n\n',APD);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% END FLAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%