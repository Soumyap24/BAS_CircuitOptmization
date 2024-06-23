
function BAS()
 clear all
close all

%antenna distance
d0=0.1;
d1=3;
d=d1;
eta_d=0.95;

%random walk
l0=0.0;
l1=0.0;
l=l1;
eta_l=0.95;

%steps
step=20;%step length
eta_step=0.95;
n=40;%iterations
k=5;%space dimension
x0=100*abs(rands(k,1));
%x0=[90;10;10;30;50;10];
x=x0;
xbest=x0;
fbest=fun(xbest);
fbest_store=fbest;
x_store=[0;x;fbest];
display(['0:','xbest=[',num2str(xbest(1)),'- ',num2str(xbest(2)),'],fbest=',num2str(fbest)])

%iteration

for i=1:n
     
    dir=rands(k,1);
    dir=dir/(eps+norm(dir));
    xleft=abs(x+dir*d);
    %xleft=enforceBiasingConstraints(xleft);
    xleft(2)=max(xleft(1)*1e3/11,xleft(2)*1e2)*1e-2;
    xleft(3)=min(xleft(3),50);
   

    fleft=fun(xleft);
    xright=abs(x-dir*d);
    %xright=enforceBiasingConstraints(xright);
    
     xright(2)=max(xright(1)*1e3/11,xright(2)*1e2)*1e-2;
     xright(3)=min(xright(3),50);

    fright=fun(xright);
    w=l*rands(k,1);
    x=abs(x-step*dir*sign(fleft-fright)+w);
     x(2)=max(x(1)*1e3/11,x(2)*1e2)*1e-2;
     x(3)=min(x(3),50);
   % x=enforceBiasingConstraints(x);
     x = max(x, 10);

    % Ensure values are within reasonable upper bounds
    upper_bound = 200; % Example upper bound for resistors
    x = min(x, 200);
  

    f=fun(x);
   
    
   
    %%%%%%%%%%%
    if f<fbest 
        xbest=x;
        fbest=f;
    end
    %%%%%%%%%%%
    x_store=cat(2,x_store,[i;x;f]);
    fbest_store=[fbest_store;fbest];
    display([num2str(i),':xbest=[',num2str(xbest(1)),' ,',num2str(xbest(2)),',',num2str(xbest(3)),',',num2str(xbest(4)),',',num2str(xbest(5)),'],fbest=',num2str(fbest)])
    %%%%%%%%%%%
    d=d*eta_d+d0;
    l=l*eta_l+l0;
    step=step*eta_step;
end



figure(3),clf(3),
% plot(x_store(1,:),x_store(4,:),'r-o')
% hold on,
plot(x_store(1,:),fbest_store,'b-.')
xlabel('iteration')
ylabel('minimum value')
end

function x = enforceBiasingConstraints(x)
    % Ensure no component value is less than 0.1 (or a suitable minimum value)
    x = max(x, 0.1);

    % Define scaling factors for the resistors
    scaling_factors = [1e3, 1e2, 1e2, 1e1, 1e4];
    
    % Scale the resistor values
    Ra = x(1) * scaling_factors(1);
    Rb = x(2) * scaling_factors(2);
    Rc = x(3) * scaling_factors(3);
    RE = x(4) * scaling_factors(4);
    R = x(5) * scaling_factors(5);

    % Define constants
    Vcc = 12; % Supply voltage
    V_BE = 0.7; % Base-emitter voltage for silicon BJT
    beta = 300; % Current gain for the transistors
    V_CE_sat = 0.3; % Saturation voltage for V_CE
    V_T = 0.026; % Thermal voltage at room temperature (26 mV)

    % Calculate base voltage V_B
    V_B = Vcc * (Rb / (Ra + Rb));
    
    % Ensure V_B is always above 0.8V
    if V_B < 1
        % Adjust Rb to ensure V_B >= 0.8V
        V_B = 1;
        Rb = V_B * Ra / (Vcc - V_B);
        x(2) = Rb / scaling_factors(2); % Update scaled Rb
    end

    % Calculate base current I_B for Q1
    I_B = (V_B - V_BE) / (RE + (Rc / beta));

    % Ensure Vcc/(Ra + Rb) is much greater than I_B
    if (Vcc / (Ra + Rb)) < 10 * I_B
        scale_factor = 10 * I_B * (Ra + Rb) / Vcc;
        Ra = Ra * scale_factor;
        Rb = Rb * scale_factor;
        x(1) = Ra / scaling_factors(1); % Update scaled Ra
        x(2) = Rb / scaling_factors(2); % Update scaled Rb
        V_B = Vcc * (Rb / (Ra + Rb));
        I_B = (V_B - V_BE) / (RE + (Rc / beta));
    end

    % Calculate collector current I_C for Q1
    I_C = beta * I_B;

%     % Ensure RE * I_C is much greater than V_T
%     if (RE * I_C) < 10 * V_T
%         RE = 10 * V_T / I_C;
%         x(4) = RE / scaling_factors(4); % Update scaled RE
%         I_B = (V_B - V_BE) / (RE + (Rc / beta));
%         I_C = beta * I_B;
%     end

    % Calculate collector-emitter voltage V_CE for Q1
    V_CE_Q1 = Vcc - I_C * Rc-RE*I_C;

    % Ensure V_CE_Q1 is above saturation voltage
    if V_CE_Q1 < V_CE_sat
        Rc = (Vcc - V_CE_sat) / I_C;
        x(3) = Rc / scaling_factors(3); % Update scaled Rc
    end

    % Calculate collector voltage V_C for Q1
    V_C_Q1 = Vcc - I_C * Rc;

    % Ensure V_B < V_C for Q1
    if V_B >= V_C_Q1
        max_Rb = Ra * (V_C_Q1 / Vcc) / (1 - V_C_Q1 / Vcc);
        Rb = min(Rb, max_Rb);
        x(2) = Rb / scaling_factors(2); % Update scaled Rb
        V_B = Vcc * (Rb / (Ra + Rb));
    end

  

    % Update the resistor values
    x(1) = Ra / scaling_factors(1);
    x(2) = Rb / scaling_factors(2);
    x(3) = Rc / scaling_factors(3);
    x(4) = RE / scaling_factors(4);
    x(5) = R / scaling_factors(5);
      % Ensure no negative or zero values for resistors again
    x = max(x, 10);

    % Ensure values are within reasonable upper bounds
    upper_bound = 200; % Example upper bound for resistors
    x = min(x, upper_bound);
end


function yout = fun(x)
    % Update netlist with new parameters
    updateNetlist(x);

    % Run LTSpice simulation
    dos('C:\Users\Soumya\OneDrive\Desktop\BAS\second_attempt\BatchCall.bat');
    pause(2); % Let LTSpice finish simulating
    dos('C:\Users\Soumya\OneDrive\Desktop\BAS\second_attempt\BatchEnd.bat'); % Closes LTSpice after the .raw file is created
    raw_data = LTspice2Matlab('example.raw');

    result = processRawData(raw_data);
     [inputSignal, outputSignal] = extractSignals(raw_data);

    % Objective function
    yout = computeObjective(result,inputSignal,outputSignal,x);
end


function yout = computeObjective(result, inputSignal, outputSignal, x)
    % Compute the output power and efficiency
    Vout_pp = result(2); % Peak-to-peak output voltage
    RL = 8; % Load resistance (ohms)
    
    % Convert peak-to-peak voltage to RMS voltage
    Vout_rms = Vout_pp / (2 * sqrt(2));
    
    % Desired output power
    Pout_desired = 1; % Desired output power in watts
    
    % Compute the actual output power
    Pout_actual = (Vout_rms^2) / RL;
    
    % Supply voltage
    Vcc = 12; % Supply voltage in volts
    
    % Compute efficiency
    %efficiency = (Pout_actual / (Vcc * I_total)) * 100; % Example efficiency calculation

    % Compute THD
    THD = computeTHD(outputSignal);
    
    % Compute linearity
    linearity = computeLinearity(inputSignal, outputSignal);
    
%     % Extract operating point values
%     V_B = Vcc * (x(2) *1e2/ (x(1)*1e3 + x(2)*1e2)); % Base voltage
%     I_B = (V_B - 0.7) / (x(4)*1e1 + (x(3)*1e2 / 300)); % Base current
%     I_C = 300 * I_B; % Collector current
%     V_C = Vcc - I_C * x(3)*1e2; % Collector-emitter voltage
%     V_CE=V_C-I_C*x(4)*1e1;
V_C=Vcc-result(3)*x(3)*1e2;
V_B=Vcc-result(4)*x(1)*1e3;

    % Penalty for saturation
    penalty = 0;
    V_CE_sat = 0.3; % Saturation voltage for V_CE
    if  V_B<1 || V_C<V_B || abs(mean(outputSignal))>0.2
        penalty = 10; % Add a penalty if in saturation or if V_B is not less than V_C
    end
    
   
   yout = 0.95* abs(1 - Vout_pp/2) +0.05*abs(1-abs(linearity)) +penalty;
    
    % Display the intermediate results for debugging
    disp(['Objective function value: ', num2str(yout)]);
end
function THD = computeTHD(signal)
    % Perform FFT on the signal to calculate THD
    L = length(signal); % Length of signal
    Y = fft(signal);
    P2 = abs(Y / L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    
    % Calculate power of fundamental and harmonics
    fundamental = P1(2); % The first element is the DC component, second is the fundamental
    harmonics = sum(P1(3:end).^2); % Sum of squares of the harmonic components
    
    % Compute THD
    THD = sqrt(harmonics) / fundamental;
end

function linearity = computeLinearity(inputSignal, outputSignal)
    % Ensure input and output signals are of the same length
    minLength = min(length(inputSignal), length(outputSignal));
    inputSignal = inputSignal(1:minLength);
    outputSignal = outputSignal(1:minLength);
   
    % Compute linearity by calculating the correlation coefficient
    correlation = corrcoef(inputSignal, outputSignal);

   
        linearity = correlation(1, 2);
end



function [inputSignal, outputSignal] = extractSignals(raw_data)
    % Extract input and output signals from the raw_data
    inputSignal = raw_data.variable_mat(raw_data.variable_name_list == "V(n004)", :); % Replace with the correct node
    outputSignal = raw_data.variable_mat(raw_data.variable_name_list == "V(vout)", :); % Replace with the correct node
end




function updateNetlist(x)
    % Define the netlist path
    netlist = ['C:\Users\Soumya\OneDrive\Desktop\BAS\second_attempt\example.net'];

    % Create the LTSpice netlist with updated parameters
%     code = ['C:\\Users\\Soumya\\OneDrive\\Desktop'...
%             '\\BAS\\second_attempt\\example.asc\r\n'...
%             'R1 Vout Vin ', num2str(x(1)), '\r\n'...
%             'C1 Vout 0 ', num2str(x(2)), '\r\n'...
%             'V1 Vin 0 SINE(0 {ampl} {freq})\r\n'...
%             '.params ampl = 1 freq = 15k\r\n'...
%             '.tran .01\r\n'...
%             '.backanno\r\n'...
%             '.end\r\n'];

% code=['C:\\Users\\Soumya\\OneDrive\\Desktop\\BAS\\second_attempt\\example.asc\r\n'...
%       'V1 N001 0 SINE(0 {ampl} {freq})\r\n'...
%       'R1 N002 N001 {R}',num2str(x(1)),'\r\n'...
%       'C1 N002 0 {C}',num2str(x(2)*1e-12),'\r\n'...
%       '.tran .1\r\n'...
%       '.params R=10k C=160p ampl=1 freq=10k\r\n'...
%       '.backanno\r\n'...
%       '.end\r\n'
%     ];
    code=['C:\\Users\\Soumya\\OneDrive\\Desktop\\BAS\\second_attempt\\example.asc\r\n'...
          'Ra N001 N005 {Ra}',num2str(x(1)*1e3),'\r\n'...
          'Rb N005 0 {Rb}',num2str(x(2)*1e2),'\r\n'...
          'Rc N001 N002 {Rc}',num2str(x(3)*1e2),'\r\n'...
         'RE N008 0 {RE}',num2str(x(4)*1e1),'\r\n'...
         'R N007 0 {R}',num2str(x(5)*1e4),'\r\n'...
         'Vcc N001 0 {Vcc}\r\n'...
         'Vin N004 0 SINE(0 {ampl} {freq})\r\n'...
         'D1 N002 N006 1N4148\r\n'...
         'D2 N006 N007 1N4148\r\n'...
         'Q1 N002 N005 N008 0 BC547B\r\n'...
       'Cin N005 N004 {Cin}\r\n'...
       'Cout Vout N003 {Cout}\r\n'...
       'RL N003 0 8\r\n'...
      'Q2 N001 N002 N003 0 2SCR293P5\r\n'...
     'Q3 0 N007 N003 0 2SAR293P5\r\n'...
      'C1 N008 0 100µ\r\n'...
      '.model D D\r\n'...
      '.lib C:\\Users\\Soumya\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
      '.model NPN NPN\r\n'...
      '.model PNP PNP\r\n'...
      '.lib C:\\Users\\Soumya\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.bjt\r\n'...
    '.model NMOS NMOS\r\n'...
    '.model PMOS PMOS\r\n'...
     '.lib C:\\Users\\Soumya\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.mos\r\n'...
     '.param Cin=10u Ra=10k Rb=1k Rc=2k RE=100 R=100k Cout=10u ampl=0.01 freq=1k Vcc=12 \r\n'...
     '.tran 0.1\r\n'...
    '.backanno\r\n'...
    '.end\r\n'

        ];


    fid = fopen(netlist, 'w+');
    fprintf(fid,code);
    fid=fclose(fid);
    % Check if file was opened successfully
    
end

% function result = processRawData(raw_data)
%     % Assuming raw_data is a structure with fields 'time', 'V(node)', and 'I(component)'
%     
%     % Extract voltage and current over the entire time series
% %     voltage = raw_data.I(C1); % replace 'node' with actual node name
% %     current = raw_data.I(R1); % replace 'component' with actual component name
% 
%      voltage1 = raw_data.variable_mat(raw_data.variable_name_list == "V(n004)", :);
%     voltage2 = raw_data.variable_mat(raw_data.variable_name_list == "V(vout)", :);
%     % Calculate peak-to-peak voltage and current
%     Vpp1 = max(voltage1) ;
%     Vpp2 = max(voltage2);
%     Opp1=-1*min(voltage1);
%     Opp2=-1*min(voltage2);
%   
%     
%     % Combine peak-to-peak voltage and current into result
%     result = [Vpp1,Vpp2,Opp1,Opp2];
%     
%     disp(result);
%     
% end

function result = processRawData(raw_data)
    % Assuming raw_data is a structure with fields 'time', 'V(node)', and 'I(component)'
    
    % Extract voltage and current over the entire time series
%     voltage = raw_data.I(C1); % replace 'node' with actual node name
%     current = raw_data.I(R1); % replace 'component' with actual component name

     voltage1 = raw_data.variable_mat(raw_data.variable_name_list == "V(n004)", :);
    voltage2 = raw_data.variable_mat(raw_data.variable_name_list == "V(vout)", :);
    I_C_raw = raw_data.variable_mat(raw_data.variable_name_list == "I(Rc)", :);
     I_B_raw = raw_data.variable_mat(raw_data.variable_name_list == "I(Ra)", :);
    % Calculate peak-to-peak voltage and current
    Vpp1 = max(voltage1) - min(voltage1);
    Vpp2 = max(voltage2) - min(voltage2);
  
    I_C=mean(I_C_raw);
    I_B=mean(I_B_raw);
    % Combine peak-to-peak voltage and current into result
    result = [Vpp1,Vpp2,I_C,I_B];
    
    disp(result);
    
end