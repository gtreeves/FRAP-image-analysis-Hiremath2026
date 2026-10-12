% script_bleachanalysis2
%
% The files processed will be the same as can be found in 
% an excel table "Bleaching_stats.xlsx".

% It runs the analysis on every file listed in bleaching_filenames.mat.

clear
close all
yesplot = false;
load bleaching_filenames % after running "script_gatherbleachfiles"


%
% Loop through each directory
%
soln = [];
for ii = 1:length(filenames)
	
	filename = filenames{ii};
	genotype = 'dl-mGFP';
    filenum = ii;
	
	data = run_analyze_bleach(filename,genotype,filenum,1,2);
	if ~isstruct(data)
		disp([data,' in:'])
		disp(filename)
		disp(' ')
	else
		soln = [soln;data];
	end
	disp(['ii = ',num2str(ii),' out of ',num2str(length(filenames))])
end



clk1 = clock;
clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2),'_',...
	num2strDU(clk1(4),2),'-',num2strDU(clk1(5),2),'-',num2strDU(round(clk1(6)),2)];
% clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2)];

% save(['Mat/',clk,'_Soln'],'soln')



