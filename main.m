% Robert Bain's final project code for COGS 268A @ UCI
% ELIF neurons with dendrite. Circuit is like an EMD / Reichardt detector.
% Includes motion detection experiments using said circuit / software model.

% Future ideas:
%   - COBA instead of CUBA (non-instantaneous current injection)
%   - Dendritic V ceiling related to synaptic equilibrium V
%   - Spine search dynamics
%   - Learning of weights
%   - Dendritic branching
%   - Active dendritic properties (eg NMDA, Ca++, Na+ spikes)
%   - Ensure bio-plausible parameters
%   - Voltage sensitive ch's that close briefly when initially ...
%     depolarized would lead to less leak AND FASTER PROPAGATION.

clear; close all; clc;

set(0, 'DefaultLineLineWidth', 2);

%% Parameters for the soma
tau = 0.010;                % membrane time constant
E_L = -0.070;               % leak potential (also resting potential)
Vth = -0.050;               % threshold potential (to produce spike)
Vreset = -0.080;            % reset potential (post-spike)
Cm = 100e-12;               % total membrane capacitance
G_L = Cm / tau;             % total membrane conductance (leak conductance)
% ELIF parameters for the soma
Vmax = 0.050;               % Clip spikes at this value of membrane potential
Delta_th = 0.005;           % Range of membrane potential for spike acceleration
% @TODO add some stochasticity
base_I = 178e-12;           % gives baseline firing rate 

%% Parameters for dendritic compartment #1
tau_d1 = 0.010;                
Cm_d1 = 100e-12;        

%% Parameters for dendritic compartment #2
tau_d2 = 0.010;                            
Cm_d2 = 100e-12;               

% resistor value coupling compartments
R_1 = 5e8;
R_2 = 5e8;

%% Simulation parameters
dt = 0.0001;                % time-step
tmax = 0.4;
t = 0:dt:tmax;               % vector of time-points
n_ts = length(t);               % number of timesteps in the simulation

% how many timesteps after feature N 'a' to fire 'b'
b_delays = -1250:50:1250;
% Get the mean FR for several different motion speeds. 
FRs = zeros(size(b_delays));
h_fig_i = 0;    % used for subplots across diff experimental delays

for h = 1:length(b_delays)
    b_delay = b_delays(h);
    
    % reset experiment specific values / arrays
    I = base_I * ones(size(t));     % vector for current at each time-point
    I_from_dends = zeros(size(t));  % total current from dend's integrated at soma
    V = E_L * ones(size(t));        % initialize the membrane potential vector for the soma
    V_d1 = E_L * ones(size(t));     % initialize the membrane potential vector for dend #1
    V_d2 = E_L * ones(size(t));     % initialize the membrane potential vector for dend #1
    G_L_d1 = Cm_d1 / tau_d1;
    G_L_d2 = Cm_d2 / tau_d2;
    spikes = zeros(size(t));        % initialize a vector to record spikes
    % less comparmental leak after sufficient depolarization
    c_leaks_1 = zeros(size(t));
    c_leaks_2 = zeros(size(t));
    c_leaks_1(1) = 1e-8; % comp_leak_1;
    c_leaks_2(1) = 1e-8; % comp_leak_2;
    
    for i = 2:n_ts             % loop through all time points

        % update dendritic compartment #2's Vm
        I_d2 = (V_d2(i-1) - V_d1(i-1)) / R_2;
        % I_d2 = I_d2 * (1 - comp_leak_2);
        % subtract current leaving d2 -> d1
        total_I_d2 = -I_d2 + G_L_d2 * (E_L - V_d2(i-1));

        if i >= 1750 + b_delay && i < 2250 + b_delay
           % feature N 'b' fires 
           % add current to dendritic compartments attached to 'b' axon
           total_I_d2 = total_I_d2 + 80e-11;
        end    
        V_d2(i) = V_d2(i-1) + dt * total_I_d2 / Cm_d2;

        % decrease outward current if sufficiently depolarized
        if V_d2(i) > 0.0 % 0 mV
            % @TODO stopped using comp_leak variables; update naming to reflect change
            % comp_leak_2 = max(0, (comp_leak_2 - 0.02));
            % @TODO make not hard-coded!
            G_L_d2 = 5e-9;
        else
            % comp_leak_2 = min(0.7, (comp_leak_2 + 0.02));
            G_L_d2 = min(1e-8, (G_L_d2 + 2e-11));
        end
        c_leaks_2(i) = G_L_d2; % comp_leak_2;

        % @TODO Fix so conduction is faster ...
        %%% Didn't work. Made the voltage level off right @ the switching point...
        % probabilistically send more current to soma conditioned on leak amount
        % if G_L_d1 == 1e-10;
        %    j = 2;
        % end
        
        I_d1 = (V_d1(i-1) - V(i-1)) / R_1;
        % @TODO no longer using this method, cleanup similar snippets too
        % I_d1 = I_d1 * (1 - comp_leak_1);
        % add d2 -> d1; subtract d1 -> soma
        total_I_d1 = -I_d1 + I_d2 + G_L_d1 * (E_L - V_d1(i-1));    

        if i >= 1750 && i < 2250
            % feature N 'a' fires
            % add current to dendritic compartments attached to 'a' axon
            % 'a' & 'b' both send axons to an inhibitory N. 
            % 'a' does disinhibition and 'b' does inhibition (ie exc & inh)
            % violation of Dale's rule or a second sign switching interneuron
            if b_delay > 0 | b_delay < -500
                total_I_d1 = total_I_d1 + 80e-11;
            else
                total_I_d1 = total_I_d1 + 80e-11 / 2;
            end
        end
        V_d1(i) = V_d1(i-1) + dt * total_I_d1 / Cm_d1;

        % decrease outward current if sufficiently depolarized
        if V_d1(i) > -0.007 % 0 mV
            % comp_leak_1 = max(0, (comp_leak_1 - 0.02));
            G_L_d1 = 5e-9;
        else
            % comp_leak_1 = min(0.7, (comp_leak_1 + 0.02));
            G_L_d1 = min(1e-8, (G_L_d1 + 2e-11));
        end
        c_leaks_1(i) = G_L_d1; % comp_leak_1;

        % next line: Forward Euler update of membrane potential
        exp_bit = Delta_th * exp((V(i-1) - Vth) / Delta_th);
        total_I = I(i) + G_L * (E_L - V(i-1) + exp_bit);
        % add current coming from dendritic compartment #1
        total_I = total_I + I_d1;
        if t == 1
            I_from_dends(i) = I_d1;
        else
            I_from_dends(i) = I_from_dends(i-1) + I_d1;
        end
        V(i) = V(i-1) + dt * total_I / Cm;
        if V(i) > Vmax              % if potential is greater than Vmax
            spikes(i) = 1;          % record the spiketime
            V(i) = Vreset;          % Reset the membrane potential
            V(i-1) = Vmax;          % add a spike on prior time-point for visual purposes
        end
    end
    
    FRs(h) = sum(spikes) / tmax;
    
    %%{
    
    if h == 1 | h == 20 | h == 28
        if h == 1
            figure()
        end
        subplot(3,3,3*h_fig_i+1)
        plot(t*1000, V*1000, 'k');
        if h == 1
            title(["(Soma)", "", "Time Vs Membrane Voltage"])
        end
        xlabel('time (ms)')
        ylabel('V_m (mV)')
        ylim([-80 80])

        subplot(3,3,3*h_fig_i+2)
        xlabel('time (ms)')
        yyaxis left
        plot(t*1000, V_d1*1000, 'k');
        ylabel('V_m (mV)')
        ylim([-80 80])
        set(gca, 'YColor','black')
        if h == 1
            title(["(Dendritic Compartment #1)", "", "Time Vs Membrane Voltage & Leak"])
        end
        yyaxis right
        ylabel('Conductance (nS)')
        set(gca, 'YColor','blue')
        ylim([40 110])
        yticks(50:10:100)
        hold on
        plot(t*1000, c_leaks_1*10e9, '--', 'Color', 'blue');
        hold off

        subplot(3,3,3*h_fig_i+3)
        xlabel('time (ms)')
        yyaxis left
        plot(t*1000, V_d2*1000, 'k');
        ylabel('V_m (mV)')
        ylim([-80 80])
        set(gca, 'YColor','black')
        if h == 1
            title(["(Dendritic Compartment #2)", "", "Time Vs Membrane Voltage & Leak"])
        end
        yyaxis right
        ylabel('Conductance (nS)')
        set(gca, 'YColor','blue')
        ylim([40 110])
        yticks(50:10:100)
        hold on
        plot(t*1000, c_leaks_2*10e9, '--', 'Color', 'blue');
        hold off
        
        h_fig_i = h_fig_i + 1;

        % figure()
        % plot(t*1000, I_from_dends, 'k');
        % title("current between soma and dendrite")
        % xlabel('time (ms)')
        % ylabel('current (pA)')

        % pause()
        % close all
    end
    
    %}
    
end

figure()
plot(b_delays/10, FRs, 'k');
title("Delay Time Vs Mean Firing Rate")
xlabel('time (ms)')
ylabel('Firing Rate (Hz)')
ylim([4 11])
yticks(5:1:10)
