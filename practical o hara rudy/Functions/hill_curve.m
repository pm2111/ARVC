function output=hill_curve(X,h, IC50) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function returning the desired hill curve for a given drug
%inputs
%  - h = hill coeff 
%  - IC50 = concentration at which channel is blocked by 50%

%output = fraction of original current still passing through
%(dimensionless)

output = 1.0 / 1+ (X/IC50)^h ;

end
