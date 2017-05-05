function output=modHH(t,X,varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hodgkin-Huxley Neuron Model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% State Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Membrane Potential V
V = X(1);
% INa gating variables
m = X(2);
h = X(3);
% IK gating variables
n = X(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional Inputs (varargin)
% Five optional inputs may be assigned:
%
% 1) flag_ode: 
%    - flag_ode = 0  -> "computed variables" output
%    - flag_ode = 1* -> dX output
%
% 2) pstim: stimulation protocol and parameters
%    - pstim = 0 -> No Stimulation  
%    - pstim = 1 -> Constant Current Stimulation
%                   - Current Amplitude (I) in uA/uF 
%                   *** input as [pstim I]
%    - pstim = 2 -> Periodic Current Stimulation 
%                   - Current Amplitude (I) in uA/uF
%                   - Stimulus Duration (td) in ms
%                   - Cycle Length (CL) in ms
%                   *** input as [pstim I duration CL] -> [2,50,0.2,100]*
%    - pstim = 3 -> Voltage-clamp
%
% When no values are provided, the default ones are used (marked with *)


% Set default values for optional inputs
optargs = {1,[2,80,0.5,100]};
newVals = cellfun(@(x) ~isempty(x), varargin);
optargs(newVals) = varargin(newVals);
[flag_ode, pstim]=optargs{:};


CL = varargin{2}(4); %cycle length
%ramp variables
t_sim = varargin{3};
number_steps = varargin{4};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other Variables
% Resting Membrane Potential
Vm = -60; %mV
% Cell capacitance
C = 1; %uF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Na+ current
% conductance
gNa = 120;  
% Nernst potential
eNa = 115; %mV
% gating variables m
alpha_m = (2.5-0.1*(V-Vm)) ./ (exp(2.5-0.1*(V-Vm)) -1);
beta_m = 4*exp(-(V-Vm)/18);
dm = alpha_m * (1-m) - beta_m * m;
% gating variables h
alpha_h = 0.07*exp(-(V-Vm)/20);
beta_h = 1./(exp(3.0-0.1*(V-Vm))+1);
dh = alpha_h * (1-h) - beta_h * h;
% 
INa = gNa * m^3 * h * ((V-Vm)-eNa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K+ current
% conductance
gK = 36;  
% Nernst potential
eK = -12; %mV
% gating variables n
alpha_n = (0.1-0.01*(V-Vm)) ./ (exp(1-0.1*(V-Vm)) -1);
beta_n = 0.125*exp(-(V-Vm)/80);
dn = alpha_n * (1-n) - beta_n * n;
% 
IK =  gK * n^4 * ((V-Vm)-eK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leak current
% conductance
gL=0.3;  
% Nernst potential
eL = 10.6; %mV
%
Ileak = gL * ((V-Vm)-eL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Procotols
switch pstim(1)
    case 0
    % 0 - No Stimulation
        Istim = 0;
        % update V -> X(1)
        dV = - (INa+IK+Ileak+Istim)/C;
    case 1
    % 0 - Constant Stimulation
        Istim = -pstim(2);
        % update V -> X(1)
        dV = - (INa+IK+Ileak+Istim)/C;
    case 2
    % 1 - Single Istim current
        amp = -pstim(2); % uA/uF      
        duration = pstim(3); % ms
        CL = pstim(4);
        trem = mod(t,CL);
        if trem <= duration; Istim = amp; else Istim = 0; end;        
        % update V -> X(1)
        dV = - (INa+IK+Ileak+Istim)/C;
    case 3
    % 2 - Voltage_Clamp
        Istim = 0;
        % update V -> X(1)
        dV = 0;
    case 4
    % Current Ramp
        Istim_max = -pstim(2);
        %number_steps = 5;
        delta_I = Istim_max/number_steps;
        I_stim_array = 0;
        for i = 1:number_steps
            I_stim_array(i+1)  = I_stim_array(end) +delta_I;
        end
        I_stim_interp = kron(I_stim_array,ones(1,round(t_sim/number_steps)+1));
        
        Istim = I_stim_interp(1,round(t+1));
        
         % update V -> X(1)
        dV = - (INa+IK+Ileak+Istim)/C;
       case 5
    % Current Ramp
        Istim_max = -pstim(2);
        %number_steps = 5;
        delta_I = Istim_max/number_steps;
        I_stim_array = Istim_max;
        for i = 1:number_steps
            I_stim_array(i+1)  = I_stim_array(end) - delta_I;
        end
        I_stim_interp = kron(I_stim_array,ones(1,round(t_sim/number_steps)+1));
        
        Istim = I_stim_interp(1,round(t+1));
        
         % update V -> X(1)
        dV = - (INa+IK+Ileak+Istim)/C;
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT
% flag = 1 -> OUTPUT is a vector with State Variables derivatives(dX)
if flag_ode==1
  dX = [dV dm dh dn]';  
  output = dX;
% flag = 0 -> OUTPUT is a vector of Computed Variables (e.g. currents):
else
  output = [Istim INa IK Ileak]';  
end
% END function
end