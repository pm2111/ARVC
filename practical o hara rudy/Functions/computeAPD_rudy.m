function APD=computeAPD_rudy(t,V,p,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to computer the AP duration at different % of repolarisation
%
% INPUTs:
% - t -> time
% - V -> membrane potential
% - p -> % of repolarisation (usually between 0 and 90)
% - flag -> if 1, it gives a figure as output, to visualise how the result 
%           is computed, e.g. start and end point
%
% The function gives as output the AP duration in ms, computed from the
% instant of max upstroke velocity to the instant in which the voltage
% goes back to a specific % of the AP amplitude.
%
% N.B. if there are multiple APs in the same trace, the fuction computes
% the duration of the first one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i100 = find(t>=100,1,'first');

% AP amplitude (max - min)
m = V(1);
[M,iM] = max(V(1:i100));


% Looking for the instant of max upstroke velocity
[~, i_start] = max( diff(V(1:iM)) ./ diff(t(1:iM)) );

%Vth = m + (1-p/100)*(M-m);
Vth = p/100 * m;

i_end = find(V(iM:end)<=Vth,1,'first') +iM; %why do we add iM? Doesnt the find() already give index at X% of repol?

APD = t(i_end) - t(i_start);

% Report
if flag==1
    figure
    plot(t,V,'k','LineWidth',2)
    hold on
    plot(t(i_start),V(i_start),'bd','MarkerSize',6,'MarkerFaceColor','b')
    plot(t(i_end),V(i_end),'mo','MarkerSize',6,'MarkerFaceColor','m')
    title('APD function report');
    legend('AP trace','start','end');
end

end
% end function