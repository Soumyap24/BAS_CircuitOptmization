function BAS()
 clear all
close all

%antenna distance
d0=0.01;
d1=3;
d=d1;
eta_d=0.99;

%random walk
l0=0;
l1=0.0;
l=l1;
eta_l=0.99;

%steps
step=1;%step length
eta_step=0.99;
n=100;%iterations
k=9;%space dimension
% x0=100*abs(rands(k,1));
x0=[7.36615;5.05005 ;2.05765 ;1.55599 ;6.47351 ;0.74808 ;33.201;26.7567;41.1706];
% x0=[18;18;18;18;18;18;18;18;3;10;50];
% x0=[37.8272 ;16.093 ;15.7512 ;15.4841 ;37.2477;14.9016;2.409;8.804 ;28.9391];
% x0=[10.1239 ,5.2225 ,2.7497 ,2.8713 ,7.0144 ,0.070036 ,33.7779,24.2158,43.053,]
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
    fleft=fun(xleft);
    xright=abs(x-dir*d);
    fright=fun(xright);
    w=l*rands(k,1);
    x=abs(x-step*dir*sign(fleft-fright)+w);
     f=fun(x);
   
    if f<fbest 
        xbest=x;
        fbest=f;
    end
   
    x_store=cat(2,x_store,[i;x;f]);
    fbest_store=[fbest_store;fbest];
    display([num2str(i),':xbest=[',num2str(xbest(1)),' ,',num2str(xbest(2)),' ,',num2str(xbest(3)),' ,',num2str(xbest(4)),' ,',num2str(xbest(5)),' ,',num2str(xbest(6)),' ,',num2str(xbest(7)),',',num2str(xbest(8)),',',num2str(xbest(9)),',','],fbest=',num2str(fbest)])
    
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

function yout = fun(x)
    % Update netlist with new parameters
    updateNetlist(x);

    % Run LTSpice simulation
    dos('C:\Users\Soumya\OneDrive\Desktop\BAS\second_attempt\BatchCall2.bat');
    pause(2); % Let LTSpice finish simulating
    dos('C:\Users\Soumya\OneDrive\Desktop\BAS\second_attempt\BatchEnd2.bat'); % Closes LTSpice after the .raw file is created

         raw_data = LTspice2Matlab('opamp2.raw');
   
    result = processRawData(raw_data);
     

    % Objective function
    yout = computeObjective(result);
end

function updateNetlist(x)
    % Define the netlist path
    netlist = ['C:\Users\Soumya\OneDrive\Desktop\BAS\second_attempt\opamp2.net'];


    code=['C:\\Users\\Soumya\\OneDrive\\Desktop\\BAS\\second_attempt\\opamp2.asc\r\n'...
'M1 N002 N002 N001 N001 CMOSP l=180n w=',num2str(x(1)*1e-6),'\r\n'...
'M2 N001 N002 N003 N001 CMOSP l=180n w=',num2str(x(1)*1e-6),'\r\n'...
'M3 N002 N006 N007 N007 CMOSN l=180n w=',num2str(x(2)*1e-6),'\r\n'...
'M4 N007 N005 N003 N007 CMOSN l=180n w=',num2str(x(2)*1e-6),'\r\n'...
'M5 N007 N008 N009 N009 CMOSN l=180n w=',num2str(x(3)*1e-6),'\r\n'...
'M6 N009 N008 N008 N009 CMOSN l=180n w=',num2str(x(4)*1e-6),'\r\n'...
'M7 N001 N003 N004 N001 CMOSP l=180n w=',num2str(x(5)*1e-6),'\r\n'...
'M8 N004 N008 N009 N009 CMOSN l=180n w=',num2str(x(6)*1e-6),'\r\n'...
'C1 N004 N003 {C1}',num2str(x(7)*1e-12),'\r\n'...
'C2 N004 0 {CL}',num2str(x(8)*1e-12),'\r\n'...
'I1 N001 N008 ',num2str(x(9)*1e-6),'\r\n'...
'V1 N001 0 2.5\r\n'...
'V2 N006 0 SINE(1.2 20m 1000 0 0 0) AC 1\r\n'...
'V3 N005 0 SINE(1.2 20m 1000 0 0 180)\r\n'...
'V4 N009 0 -2.5\r\n'...
'.model NMOS NMOS\r\n'...
'.model PMOS PMOS\r\n'...
'.lib C:\\Users\\Soumya\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.mos\r\n'...
'.ac dec 1000 1k 10000k\r\n'...
'.include tsmc018.lib\r\n'...
';tran 0 0.1 0.05\r\n'...
'.param   I=50u C1=3p CL=10p \r\n'...
'.backanno\r\n'...
'.end\r\n'

        ];
% * C:\Users\Soumya\OneDrive\Desktop\BAS\second_attempt\opamp2.asc
% M1 N002 N002 N001 N001 CMOSP l=180n w={L1}
% M2 N001 N002 N003 N001 CMOSP l=180n w={L2}
% M3 N002 N006 N007 N007 CMOSN l=180n w={L3}
% M4 N007 N005 N003 N007 CMOSN l=180n w={L4}
% M5 N007 N008 N009 N009 CMOSN l=180n w={L5}
% M6 N009 N008 N008 N009 CMOSN l=180n w={L6}
% I1 N001 N008 {I}
% V1 N001 0 2.5
% V2 N006 0 SINE(1.2 20m 1000 0 0 0) AC 1
% V3 N005 0 SINE(1.2 20m 1000 0 0 180)
% V4 N009 0 -2.5
% M7 N001 N003 N004 N001 CMOSP l=180n w={L7}
% C1 N004 N003 {C1}
% C2 N004 0 {CL}
% M8 N004 N008 N009 N009 CMOSN l=180n w={L8}
% .model NMOS NMOS
% .model PMOS PMOS
% .lib C:\Users\Soumya\OneDrive\Documents\LTspiceXVII\lib\cmp\standard.mos
% .ac dec 1000 1k 10000k
% .include tsmc018.lib
% ;tran 0 0.1 0.01
% .param L1=1600n L2=1600n L3=1600n L4=1600n L5=1600n L6=1600n L7=1600n L8=1600n I=50u C1=3p CL=10p
% .backanno
% .end



    fid = fopen(netlist, 'w+');
    fprintf(fid,code);
    fid=fclose(fid);
    % Check if file was opened successfully
end

function result = processRawData(raw_data)
    % Extract frequency and AC analysis data
    freq = raw_data.freq_vect;
    gain = abs(raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n004)"), :) ./ raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n006)"), :));

    % Ensure NMOS and PMOS transistors are in saturation
    Vgs_m3 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n006)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n007)"), :);
    Vds_m3 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n002)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n007)"), :);
    Vgs_m4 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n005)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n007)"), :);
    Vds_m4 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n003)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n007)"), :);
    Vgs_m5 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n008)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n009)"), :);
    Vds_m5 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n007)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n009)"), :);
    Vgs_m6 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n008)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n009)"), :);
    Vds_m6 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n008)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n009)"), :);
    Vgs_m8 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n008)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n009)"), :);
    Vds_m8 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n004)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n009)"), :);

    Vsg_m1 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n001)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n002)"), :);
    Vsd_m1 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n001)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n002)"), :);
    Vsg_m2 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n001)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n002)"), :);
    Vsd_m2 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n001)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n003)"), :);
    Vsg_m7 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n001)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n003)"), :);
    Vsd_m7 = raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n001)"), :) - raw_data.variable_mat(strcmp(raw_data.variable_name_list, "V(n004)"), :);

    Vth_NMOS = 0.7; % Assuming a threshold voltage of 0.7V for all NMOS transistors
    Vth_PMOS = -0.7; % Assuming a threshold voltage of -0.7V for all PMOS transistors

    nmos_saturation_check = all([Vds_m3 > (Vgs_m3 - Vth_NMOS), Vds_m4 > (Vgs_m4 - Vth_NMOS), Vds_m5 > (Vgs_m5 - Vth_NMOS), Vds_m6 > (Vgs_m6 - Vth_NMOS),Vds_m8 > (Vgs_m8 - Vth_NMOS)]);
    pmos_saturation_check = all([Vsd_m1 > (Vsg_m1 - abs(Vth_PMOS)), Vsd_m2 > abs(Vsg_m2 - abs(Vth_PMOS)),Vsd_m7 > abs(Vsg_m7 - abs(Vth_PMOS))]);

    disp(nmos_saturation_check)
    disp(pmos_saturation_check)

    penalty = 0;
    if nmos_saturation_check && pmos_saturation_check
        penalty = 1e6; % Apply a large penalty if transistors are not in saturation
    else penalty=0;
    end

    % Calculate DC gain (gain at the lowest frequency)
    dc_gain = gain(1);

    % Calculate bandwidth (frequency at which gain drops by 3 dB from DC gain)
    target_gain = dc_gain / sqrt(2);
    if isempty(find(gain <= target_gain, 1))
        bandwidth = 0;
        penalty = penalty + 1e6; % Additional penalty if bandwidth calculation fails
    else
        idx = find(gain <= target_gain, 1);
        bandwidth = freq(idx);
    end

    % Ensure result has three elements: [dc_gain, bandwidth, penalty]
    result = [dc_gain, bandwidth, penalty];
   disp(['DCgain: ', num2str(dc_gain)]);
    disp(['bandwidth: ', num2str(bandwidth)]);
    disp(['peanlaty ', num2str(penalty)]);
end


function yout = computeObjective(result)
    % Define the objective function considering DC gain and bandwidth
    dc_gain = result(1);
    bandwidth = result(2);
    penalty = result(3);

    % Objective function to maximize gain and bandwidth
    yout = -(dc_gain) + penalty; % Assuming we want to maximize both gain and bandwidth
    disp(['Objective function value: ', num2str(yout)]);
    disp(['DCgain: ', num2str(dc_gain)]);
    disp(['bandwidth: ', num2str((bandwidth*dc_gain)*1e-6)]);
    disp(['peanlaty ', num2str(penalty)]);
    
end