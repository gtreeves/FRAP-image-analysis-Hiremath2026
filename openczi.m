function y = ftn_FRAP(P,t,ycyt,experimenttype)
%
%
%
%
% 

c0 = P(1);
k_in = P(2);
k_out = P(3);

if ~exist('experimenttype','var')
	experimenttype = 'bleach';
end

if strcmp(experimenttype,'photoactivation')
	error('have not derived the equation for this yet')
else
	yI = ycyt.*exp(k_out*t);
	I = cumtrapz(t,yI);
	
	y = exp(-k_out*t).*(c0 + k_in*I);
end