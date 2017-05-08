% HDF5 reader for 1D fibre sims 

%startpath= '/home/scratch/testoutput/OUT_sim_duration_9000.00_BCL_600.00/' ;
startpath= '/home/scratch/OUT_BCL300.00/' ;

endpath = '/results.h5';
endpath1 = 'AbstractCardiacProblem_mSolution.h5';
subdirs = dir(startpath);
b =size(subdirs);
a = b(1);
results = zeros(a(1)-3,4);
activation_map = cell(1);
for i = 1:(a-3)
   
    midpath = subdirs(3+i,1).name;
    
    file_id = H5F.open(strcat(startpath , midpath , endpath)) ;

    
    k1_frac = str2double(midpath(end-3:end));
    na_frac = str2double(midpath(end-8:end-5));

    %data = hdf5read(filename);

    hinfo = hdf5info(strcat(startpath , midpath , endpath));

    %sim_time = hdf5read(hinfo.GroupHierarchy.Datasets(4));
    sim_time = hdf5read(hinfo.GroupHierarchy.Datasets(2));

    %activation_time = hdf5read(hinfo.GroupHierarchy.Datasets(5));
    activation_time = hdf5read(hinfo.GroupHierarchy.Datasets(3));

    %vol_phi = hdf5read(hinfo.GroupHierarchy.Datasets(1));
    vol_phi = hdf5read(hinfo.GroupHierarchy.Datasets(1));

    voltage = vol_phi(1,:,:);
    voltage_field_pace{i} =voltage(1,1,:);
    
    voltage_field_myo{i} =voltage(1,30,:);
    %the h5 file does not contain the cell coords
    
    cell_position = linspace(0,2,51);

    x_data = cell_position;
    nb = size(activation_time);
    
    voltage_activation = activation_counter(sim_time, voltage,300.00,0.2);
    interval = sim_time(end)/nb(3);
    for j =1:nb(3)
        
        %y_data = activation_time(1,:,j);
        y_data = voltage_activation(j,:);
        
        [xfit, yfit, vel_arvc, err_arvc] = fit_line(x_data(6:end-1),y_data(6:end));  %m is the vel in cm/ms ( *10 for m/s))

        %xfit, yfit, vel_arvc, err_arvc] = fit_line(x_data(15:end),y_data(15:end));  %m is the vel in cm/ms ( *10 for m/s))
        %[vel_arvc] = fit_rough_line(x_data(5:end),y_data(5:end)); %use if first method fails
        vel_arvc = 10 * vel_arvc; %in m/s
        %[xfit1, yfit1, vel_control, err_control ] = fit_line(x_data1,y_data1);  %m is the vel in cm/ms ( *10 for m/s))
        if abs(vel_arvc)< 0.05
            vel_arvc = 0.0;
        end
        if abs(vel_arvc)> 1.0
            vel_arvc = 0.0;
        end
        if vel_arvc < 0 
            vel_arvc =0.0;
        end
        vel_field{i,j} = vel_arvc;
        activation_map{i,j} = y_data;
        
        min_voltage(j) = mean(min(voltage(1,6:end,1+(j-1)*interval:(j)*interval),[],3));
        
    end
    %classify resting potential 
     mean_min_pot_cell = mean(min_voltage);
    %find mean velocity of fibre over many beats
     a = vel_field{i,5:end};
     avg=sum(a(a~=0))/sum(a~=0);

   
    results(i,1) = na_frac;
    results(i,2) = k1_frac;
    results(i,3) =  vel_arvc;
    results(i,4) = mean_min_pot_cell;

    vel_matrix(na_frac*10, k1_frac*10) = avg;
    min_pot_matrix(na_frac*10, k1_frac*10) = mean_min_pot_cell;
    activation(i,:) = y_data;
    %%%%%%%%find gradient of curve, excluding initial cells 

end






%Show conduction block at low scaling factors

figure()
hold all
scatter(sim_time,voltage_field_pace{1,85},'DisplayName','pacemaker cell')
scatter(sim_time,voltage_field_myo{1,85},'DisplayName','mid-fibre myocyte')
set(gca,'FontSize',25,'fontWeight','bold')
legend('show')
xlabel('time [ms]')
ylabel('membrane potential [mV]')





figure()
hold all;

for i = 1:(a-4)
plot( cell_position,activation_map{i+1},'bl','LineWidth', 2.0, 'DisplayName', 'ARVC')
%plot( activation_ARVC(:,2),activation_ARVC(:,1),'o','LineWidth', 2.0, 'DisplayName', 'ARVC Na and K1 complex')
%plot( xfit,yfit, 'LineWidth', 2.0, 'DisplayName', 'fit');
%plot(xfit1,yfit1, 'LineWidth', 2.0, 'DisplayName', 'fit2');
xlabel('cell distance from excitation edge (cm)');
ylabel('activation time (ms) ');
%legend('show');
set(gca,'FontSize',22,'fontWeight','bold')

end



%plot activation map in 2D
figure()
hold all;


imagesc( activation(30,:))
%plot( activation_ARVC(:,2),activation_ARVC(:,1),'o','LineWidth', 2.0, 'DisplayName', 'ARVC Na and K1 complex')
%plot( xfit,yfit, 'LineWidth', 2.0, 'DisplayName', 'fit');
%plot(xfit1,yfit1, 'LineWidth', 2.0, 'DisplayName', 'fit2');
xlabel('cell distance from excitation edge (cm)');
c = colorbar();
%legend('show');
set(gca,'FontSize',30,'fontWeight','bold')
labels = [];
set(gca, 'yTickLabel', labels); % Change x-axis ticks labels.
xlim([1 51]);
ylabel(c,'activation time [ms]','Fontsize',30,'fontWeight', 'bold');



figure()

%plot2( cell_position,activation_map{i},'LineWidth', 2.0, 'DisplayName', 'ARVC')
%surf(results(i,1), results(i,2), ve)
imagesc(vel_matrix);
c= colorbar;
%plot( activation_ARVC(:,2),activation_ARVC(:,1),'o','LineWidth', 2.0, 'DisplayName', 'ARVC Na and K1 complex')
%plot( xfit,yfit, 'LineWidth', 2.0, 'DisplayName', 'fit');
%plot(xfit1,yfit1, 'LineWidth', 2.0, 'DisplayName', 'fit2');
xlabel('Scaling factor of GK1');
ylabel('Scaling factor of GNa');
%ylabel(c,'conduction velocity [m/s]','Fontsize',25,'fontWeight', 'bold');
%legend('show');
title('CV [m/s] measurement for BCL 600ms and healthy [K+]ext concen');
set(gca,'FontSize',30,'fontWeight','bold')
labels = [ 0.2  0.4 0.6  0.8  1.0];
set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
set(gca, 'yTickLabel', labels); % Change x-axis ticks labels.




figure()

%plot2( cell_position,activation_map{i},'LineWidth', 2.0, 'DisplayName', 'ARVC')
%surf(results(i,1), results(i,2), ve)
imagesc(min_pot_matrix);
c= colorbar;
%plot( activation_ARVC(:,2),activation_ARVC(:,1),'o','LineWidth', 2.0, 'DisplayName', 'ARVC Na and K1 complex')
%plot( xfit,yfit, 'LineWidth', 2.0, 'DisplayName', 'fit');
%plot(xfit1,yfit1, 'LineWidth', 2.0, 'DisplayName', 'fit2');
xlabel('Fraction of IK1 channel open');
ylabel('Fraction of INa channel open');
ylabel(c,'mean resting potential of the fibre [mV]','Fontsize',22,'fontWeight', 'bold');
%legend('show');
set(gca,'FontSize',22,'fontWeight','bold')
labels = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
set(gca, 'yTickLabel', labels); % Change x-axis ticks labels.












colororder = [
	0.00  0.00  1.00
	0.00  0.50  0.00 
	1.00  0.00  0.00 
	0.00  0.75  0.75
	0.75  0.00  0.75
	0.75  0.75  0.00 
	0.25  0.25  0.25
	0.75  0.25  0.25
	0.95  0.95  0.00 
	0.25  0.25  0.75
	0.75  0.75  0.75
	0.00  1.00  0.00 
	0.76  0.57  0.17
	0.54  0.63  0.22
	0.34  0.57  0.92
	1.00  0.10  0.60
	0.88  0.75  0.73
	0.10  0.49  0.47
	0.66  0.34  0.65
	0.99  0.41  0.23
];
%ROLE OF IK1
figure()
set(gca, 'ColorOrder', colororder);


hold all;
for i =1:10
    plot(linspace(0,1,10),vel_matrix(i,:),'-o','LineWidth', 2.0,'DisplayName',strcat('Scaling factor of GNa = ',num2str(results(i*10-1,1))));
end
legend('show');
xlabel('Scaling factor of GK1');
ylabel('Conduction velocity [m/s]');
set(gca,'FontSize',30,'fontWeight','bold')


%ROLE OF INa
figure()
set(gca, 'ColorOrder', colororder);


hold all;
for i = [2 5 9]
    plot(linspace(0,1,10),vel_matrix(:,i),'-o','LineWidth', 2.0,'DisplayName',strcat('Scaling factor of GK1 = ',num2str(results(i,2))));
end
legend('show');
xlabel('Scaling factor of GNa');
ylabel('Conduction velocity [m/s]');
set(gca,'FontSize',30,'fontWeight','bold')

%Resting Potential and IK1
%activation_counter(sim_time, voltage,300, 0.2);




figure()
hold all

for i=1:50
scatter(sim_time(1:3000),voltage(1,i,1:3000))
end





function[times ] = activation_counter(sim_time, voltage,period, step_size)
%%% NB> this script is dependent on the timestep of the simulation
%period = 600.0;
%step_size= 0.2;
cell_nr = size(voltage);
cnr = cell_nr(2);
beat_nr = cell_nr(3)*step_size/period;
%stim_time_ind = find( abs(voltage(1,5,1:1500) - 40) < 0.5);
 
interval = (period/step_size );

reject_time = 5;%reject if time taken is larger that 10 ms between adjacent nodes
reject_index = reject_time/step_size;
% times = zeros(beat_nr,cnr);
% peaks = zeros(beat_nr, cnr);
% indices = zeros(beat_nr,cnr);
for j=1:beat_nr
    %[stim_max, stim_time_ind] = max( voltage(1,5,x:y));
    x = 1+ (j-1)*interval;
    y = x+ reject_index;
    for i=1:cnr-1
   
       %ind = find( abs(voltage(1,i,x:y)-25) < 20);
       [stim_max, ind] = max( diff(voltage(1,i,x:y)));
       ind = x +ind; %we start with a small sample of the total array
      % if stim_max < 50 && stim_max >5.0
           times(j,i) = sim_time(ind);
            peaks(j,i) = stim_max;
            indices(j,i) = ind;
     %  else
%            times(i)=0;
%            peaks(i) =stim_max;
%            indices(i) = ind;
           %break
      % end
       %ind = ind(1);
      %  if (abs(argvalue) - abs(stim_max)) > 20
        %    times(j,i) = 0;
        %    break
      %  else
            
       % end
        x = ind-1;
        y = x + reject_index;
        disp(x);
        disp(y)
    end
end

end




function [xf , yf, vel, err] =  fit_line(x,y)
coeffs = polyfit(x, y, 1);
%get fitted values
fittedX= linspace(min(x),max(x),200);
fittedY = polyval(coeffs, fittedX);
xf = fittedX;
yf = fittedY;
grad = coeffs(1);
vel = 1.0/grad; 
a =size(xf) ;
var = 1/a(2) * sum((y - polyval(coeffs, x)).^2 ) ;  
err = vel * sqrt( var/grad^2);
end


    
function [    vel ] = fit_rough_line(x,y)
%NB very crude approx, cannot guarantee accuracy
z = y(end) - y(1);
g = z/x(end);
vel = 1/g;
end

