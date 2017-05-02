% HDF5 reader for 1D fibre sims 
close all
clear all
startpath= '1D/' ;
endpath = '/results.h5';

subdirs = dir(startpath);
b =size(subdirs);
a = b(1);
results = zeros(a(1)-3,4);
activation_map = cell(1);
vel_field = cell(10,10);
voltage_field = cell(10,10);
for i = 1:(a-3)
    disp(strcat(i, '% complete'));
    midpath = subdirs(3+i,1).name;
    
    file_id = H5F.open(strcat(startpath ,midpath, endpath)) ;
     
    
    k1_frac = str2double(midpath(end-3:end));
    na_frac = str2double(midpath(end-8:end-5));

    %data = hdf5read(filename);

    hinfo = hdf5info(strcat(startpath  , midpath,endpath));

    sim_time = hdf5read(hinfo.GroupHierarchy.Datasets(2));

    activation_time = hdf5read(hinfo.GroupHierarchy.Datasets(3));
    
    vol_phi = hdf5read(hinfo.GroupHierarchy.Datasets(1));
    voltage = vol_phi(1,:,:);
    %the h5 file does not contain the cell coords

    cell_position = linspace(0,2,51);

    x_data = cell_position;
    l = size(activation_time);
    for k = 1:l(3)
    y_data = activation_time(1,:,k);

    [xfit, yfit, vel_arvc, err_arvc] = fit_line(x_data(15:end),y_data(15:end));  %m is the vel in cm/ms ( *10 for m/s))
    %[vel_arvc] = fit_rough_line(x_data(5:end),y_data(5:end)); %use if first method fails
    if abs(vel_arvc) > 100
        vel_arvc = 0.0;
    end
    vel_mat(k) = 10 * vel_arvc; %in m/s
    
    %[xfit1, yfit1, vel_control, err_control ] = fit_line(x_data1,y_data1);  %m is the vel in cm/ms ( *10 for m/s))
    
    end
    %classify resting potential 
    vel_field{na_frac*10,k1_frac*10} = vel_mat;
    voltage_field{na_frac*10,k1_frac*10} = voltage;
     %mean_min_pot_cell = mean(min(voltage(1,7:end,:),[],3));
    %mean_min_pot_cell = mean(min(voltage(1,7:end,end-2000:end),[],3));
    mean_min_pot_cell = mean(min(voltage(1,7:end,:),[],3));

    activation_map{end+1} = [y_data];
    results(i,1) = na_frac;
    results(i,2) = k1_frac;
    results(i,3) =  vel_arvc;
    results(i,4) = mean_min_pot_cell;

    vel_matrix(na_frac*10, k1_frac*10) = vel_arvc;
    min_pot_matrix(na_frac*10, k1_frac*10) = mean_min_pot_cell;
    activation(i,:) = y_data;
    %%%%%%%%find gradient of curve, excluding initial cells 

end


%compute avg vel for each channel fraction combination
for i = 1:10
    for j=1:10
        z = vel_field{i,j};
        if size(vel_field{i,j}) >0
            
            z = z(3:end); %only take into account the from 10th beat on
            %CHANGE THE ABOVE FOR 30 beat sim to
           % z = z(10:20); %only take into account the from 10th beat on

            avg=sum(z(z~=0))/sum(z~=0);
            vel_mat_ave(i,j) = avg;
            dummy = find(vel_field{i,j} ==0);
            %nr_beats(i,j) = dummy(1) -1;

        end
        

    end
end

figure()
hold all
scatter(sim_time,voltage(1,1,:));
scatter(sim_time,voltage(1,30,:));



figure()
imagesc(vel_mat_ave);
c = colorbar();
%legend('show');
xlabel('Scaling factor of GK1');
ylabel('Scaling factor of GNa');
%ylabel(c,'conduction velocity [m/s]','Fontsize',30,'fontWeight', 'bold');
%legend('show');
set(gca,'FontSize',27,'fontWeight','bold')
labels = [ 0.2  0.4  0.6  0.8 1.0];
set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
set(gca, 'yTickLabel', labels); % Change x-axis ticks labels.
title('CV [m/s] measurement for BCL 300ms')





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


imagesc( activation(1,:))
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
xlabel('Fraction of IK1 channel open');
ylabel('Fraction of INa channel open');
ylabel(c,'conduction velocity [m/s]','Fontsize',25,'fontWeight', 'bold');
%legend('show');
set(gca,'FontSize',25,'fontWeight','bold')
labels = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
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
for i =3:7
    plot(linspace(0,1,10),vel_mat_ave(i,:),'-o','LineWidth', 2.0,'DisplayName',strcat('INa fraction =',num2str(results(i*10-1,1))));
end
legend('show');
xlabel('Fraction of IK1 channel open');
ylabel('Conduction velocity [m/s]');
set(gca,'FontSize',30,'fontWeight','bold')


%ROLE OF INa
figure()
set(gca, 'ColorOrder', colororder);


hold all;
for i = [1 5 10]
    plot(linspace(0,1,10),vel_mat_ave(:,i),'-o','LineWidth', 2.0,'DisplayName',strcat('IK1 fraction =',num2str(results(i,2))));
end
legend('show');
xlabel('Fraction of INa channel open');
ylabel('Conduction velocity [m/s]');
set(gca,'FontSize',30,'fontWeight','bold')

%Resting Potential and IK1






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

