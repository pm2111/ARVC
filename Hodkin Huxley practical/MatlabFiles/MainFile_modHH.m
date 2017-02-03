close all; clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main File to run the Hodgkin-Huxley Neuron Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials = 100;
Y =cell(trials);
res1 = cell(1);
res2 = cell(1);

res1{1}= ['Amplitude'];
res2{1} = ['Cyclic?'];

for i = 1:trials

%% Main Settings: 
% model_name: 
mod = @modHH;
% ODE settings
t_sim = 500.; % ms
ODEstep = 0.1; % ms 
options=odeset('MaxStep',ODEstep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Conditions,e
% [V_0, m_0, h_0, n_0];
CI = [-60, 0.5, 0.5, 0.5];
disp "works"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other inputs
input_args ={};
% flag_ode: 1 when solving ODEs, 0 when computing variables
input_args{1} = 1;
% modify input args
stim_amplitude = 270 -80*i/30;
input_args{2} = [2,stim_amplitude,t_sim,100];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulations
[t,y] = ode15s(mod,[0 t_sim],CI,options,input_args{:});
 %save the final conditions of the current sim and re start with these
CI = y(end,:);
Y{i} = y;
% Computed Variables
input_args{1}=0; 
lCVs=size(feval(mod,t(1),y(1,:),input_args{:}),1);
CVs = zeros(length(t),lCVs);
for j=1:size(y,1)
    CVs(j,:)=feval(mod,t(j),y(j,:),input_args{:});
end


if( Y{i}(end-1,1) - Y{i}(end-10,1) < 0.0001)
   res1{end+1} = stim_amplitude; 
   res2{end+1} =0 ;

else
   res1{end+1} = stim_amplitude; 
   res2{end+1} =1 ;

  
end
end
% fileID = fopen('amplitude of stim.txt','w');
% fprintf(fileID,'%6s  %12s\n','I (microA)','state');
% fprintf(fileID,'%6.2f %12.8f\n',res);
% 
% fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures:
% Figure 1: Membrane Potential
cc=hsv(12); %create a diverse colour palette

h1=figure(); hold on;

for k =1:10
    k = k*10;
    plot(t(1:1000),Y{k}(1:1000,1),'LineWidth',1.5, 'Color',cc(k/10,:) , 'DisplayName', strcat('Stimulus amplitude (micro A)', num2str(270 - 80*k/30))); 

    %[the_max, index_of_max] = max(Y{k}(1:500,:,1));
    %index_of_max(1)

end
ylabel('Membrane Potential (mV)'); xlabel('time (ms)');
%set(legend_handle,'Interpreter','latex')
legend('show')
grid on


%Figure 2: Gating Variables
h2=figure(2); 
subplot(3,1,1);
plot(t,y(:,2),'LineWidth',2); 
ylabel('gate m'); xlabel('time (ms)');
subplot(3,1,2);
plot(t,y(:,3),'LineWidth',2); 
ylabel('gate n'); xlabel('time (ms)');
subplot(3,1,3);
plot(t,y(:,4),'LineWidth',2); 
ylabel('gate h'); xlabel('time (ms)');
%Figure 3: Ionic Currents
h3=figure(3); 
subplot(2,2,1);
plot(t,CVs(:,1),'LineWidth',2); 
ylabel('I_{stim} (uA/uF)'); xlabel('time (ms)');

subplot(2,2,2);
plot(t,CVs(:,2),'LineWidth',2); 
ylabel('I_{Na} (uA/uF)'); xlabel('time (ms)');

subplot(2,2,3);
plot(t,CVs(:,3),'LineWidth',2); 
ylabel('I_{K} (uA/uF)'); xlabel('time (ms)');

subplot(2,2,4);
plot(t,CVs(:,4),'LineWidth',2); 
ylabel('I_{leak} (uA/uF)'); xlabel('time (ms)');
%% END FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%