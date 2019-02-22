%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First attempt toward a script for controlling 
% the testbench for the STEP-model
% 
% Initial Authors: Johan Raman, Pieter Rombouts
% Author of the current version: Annelien Dewagtere
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisations, data of the assignment etc.

  close all
  clear all

  MyName = 'Dewagtere_Annelien'; 
      % use underscores, given name at the end
  Experiment = 'STEP'  % valid choices: 'VCO' 'PLL' 'XY' 'STEP' 'FSK'
  modelName = strcat(Experiment,'_TB'); % build modelname
  
  f0_desired = 100*10^3
    % desired free running frequency
  Q_desired = 1/sqrt(2)
    % desired closed-loop Q-factor
  VCC = 10
    % Supply Voltage [V]

% End of Initialisations, data of the assignment etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal values of the LM656.
% These values are taken from the schematic and you can use them.
% All resistances are in [Ohm]

  R1  = 7.2e3; R2  = R1; R4  = 200; R8  = 8.1e3; R10 = 1.75e3;
  R11 = 3.8e3; R12 = 3.6e3; R14 = 1e3; R15 = 205; R19  = 6.5e3;
  R20  = 4.7e3; R22  = 4.3e3; R25  = 2.6e3; R26  = 8.4e3;

% END OF internal values of the LM656.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived values of the LM656 - VCO. 
% It is possible to calculate these values from the schematic.
%         Part 1: given quantities

  %%%%%%%%%%%%%%%%%%%%%%
  % Some biasing values
  %%%%%%%%%%%%%%%%%%%%%%
  VC = 2*VCC;  % total supply voltage +VCC -(-VCC) = 2 VCC
  VD = 0.6;    % diode trheshold Voltage [V]
  Ibias  = (VC-VD)/(R4+R8+R10+R11)  % bias current [A]

  %%%%%%%%%%%%%%%%%%%%%%
  % Schmitt trigger
  %%%%%%%%%%%%%%%%%%%%%%
  Vthlow = (R25-R19)/(R25+R19)*VCC + (R25+3*R19)/(R25+R19)*VD;
     % low threshold voltage of the Schmitt trigger = [gift for you]
  Vhigh  = VCC-VD
     % High output voltage of the Schmitt-trigger.
     % This is also the high level of the square wave that is the VCO output
     % signal.
  Vlow   = (1/R19+1/R20-1/R25)/(1/R19+1/R20+1/R25)*VCC+(1/R26+1/R25-1/R19)/(1/R19+1/R20+1/R25)*VD-VD
     % Low output voltage of the Schmitt-trigger.
     % This is also the lowlevel of the square wave that is the VCO output
     % signal.
  Vthhigh = Vlow+3*VD
     % high threshold voltage of the Schmitt trigger = [gift for you]
  Vswitch = VCC-R22*Ibias    
  % reference-voltage to which Vsquare is compared to switch the current direction
    
% END OF Derived values of the LM656 - VCO. -- Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived values of the LM656 - VCO. 
% It is possible to calculate these values from the schematic.
%         Part 2: Values for you to calculate 

  V7bias = VCC-Ibias*R12/2 
    % the quiescent voltage at pin 7 [Volt]
  Ro7 = R12
    % the output resistance at pin 7 [Ohm]
  Ro6 = 1 / (1 / R10 + 1 / (R11 + R8 + R4))
    % the output resistance at pin 6 [Ohm]
    
 % ADD MORE IF NEEDED.

% END OF Derived values of the LM656 - VCO. Part 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% external values LM656 - VCO.
% In normal operation these components values must be sized
% by the user.
% All resistances are in [Ohm]

  R_T =464
    % "timing" resistor to set f0 [Ohm]
  C_T = 6.78 *10^-9
  %RoundESeries(R12*24/(10*2*(R10+R11+R8+R4))/f0_desired/R_T,6)
    % "timing" capacitor to set f0 [Farad]
    % optionally: use RoundESeries(C_T,6) to round to and E6 series value
    
  f0 = R12*24/(10*2*(R10+R11+R8+R4))/(R_T*C_T)
     % In a good design, obviously this should be equal to f0_desired
  KV = -pi*f0*24/(5*VC*R12*24/(10*2*(R10+R11+R8+R4)))
    % in [rad/s /V]
    
  RKp = RoundESeries((Ro6-(Ro6+Ro7)*2*f0*pi*R14*2/(KV*(2/Q_desired^2/pi+1/2)*-2*R12*VD))/(2*f0*pi*R14*2/(KV*(2/Q_desired^2/pi+1/2)*-2*R12*VD)-1), 12)
 % RKp=  RoundESeries(10^24,12)
    % optional Kp tuning resistor [Ohm]
    % optionally: use RoundESeries(Rkp,12) to round to and E12 series value
  RLP = 1 / (1 / Ro7 + 1 / (Ro6 + RKp))
    % effective loop filter resistance [Ohm]
    % This is actually a derived parameter (and therefore it belongs it 
    % in the next section), but you may already
    % need it for the sizing of CLP

% END OF external values LM656 - VCO.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived VCO parameters
% These are the external VCO-parameters

  MaxDeltaV7 = (R12 / R14) * (VD / 2)*(Ro6+RKp)/(Ro6+Ro7+RKp)
    % [V] maximal voltage swing on pin 7
  KP =-2*MaxDeltaV7/pi
    % [V/rad]
  fH   = KV * KP / 4
     % [Hz] Hold range
  CLP = RoundESeries(Q_desired^2 /(RLP*KV*KP ), 6);
    % capacitance used in the lowpass loop filter [F]
    % optionally: use RoundESeries(CLP,6) to round to and E6 series value
  tauLP = RLP*CLP
    % time constant of the lowpass loop filter [s]
  fLP =  1/(tauLP*2*pi)         
    % [Hz] 
    % cut off frequency of the lowpass loop filter  
  fC   = sqrt(fH * fLP)
     % [Hz] capture range

 % ADD MORE IF NEEDED.

% ENd of Derived VCO parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the testbench parameters & simulation
% + running the simulation + generating the output

    % STEP test bench parameters & simulation
    fin =100    % [Hz]
       % frequency of input square wave (used as approximaton of a step) 
    A   = 0.5
       % input step amplitude in Volt
    DelayTime = 1e-4
       % delay time [s]
       % set to something large enough to ensure that the start-up transient of
       % your PLL is finished
    max_frequentie = 2*f0;
    TimeStep = 1/10/max_frequentie;
    SimulationTime = 3e-4  % [s]
       % set to something of the order of the DelayTime + 10 times the expected rise time

    tic  % start built in matlab timer
    tmp = sim(modelName,'StopTime', num2str(SimulationTime), ...
              'MaxStep',num2str(TimeStep), 'ZeroCross','on');
    sig_in = tmp.get('sig_in'); sig_out = tmp.get('sig_out'); 
    sig_out2 = tmp.get('sig_out2'); clear tmp;
    elapsedTime = toc  % stop built in matlab timer
    
    %generate some warnings
    if (elapsedTime > 600), % if the simulation takes more than 10 minutes
      disp('Do not even think of submitting this script. Your professor will kill you!');
      disp('Modidify your settings such that the simulation time is less long.');
    elseif (elapsedTime > 120), % if the simulation takes more than 2 minutes
      disp('If you run this simulation on a fast computer, then probably you should reconfigure some settings');
      disp('prior to submitting your files. Such a long simulation time is not acceptable for the submitted version.');
      disp('Try to make your simulation run faster.')
    end
    
    %generate output
    figure(1)
    plot(sig_out.time,sig_out.signals.values)
    hold on
    plot(sig_in.time, sig_in.signals.values,'r')
    plot(sig_out2.time,sig_out2.signals.values,'k','LineWidth', 3)


% ENd of Setting up the PLL testbench etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate correct files to submit to Minerva
% 
  if 0, % if you don't want to create the zip file: alter to 'if 0,'
    GenerateZipFile(MyName,Experiment);
  end % if
    
% ENd of Generate correct files to submit to Minerva
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
