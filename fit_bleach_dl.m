function data = analyze_bleach_czi(filename,genotype,filenum,data_ch,varargin)
%Extracts nuc & cyt time course intensities; segments and tracks nuclei and
%surrounding cytoplasm.

%function data = analyze_bleach(filename,yesplot)
%
% "filename": filename (incl path) to the lsm file that contains the
%	bleaching time course.
% "genotype": string containing the genotype of the embryo in question
% "data_ch": the numerical tag for the channel where the bleaching data are
%
% Optional argument varargin can consist of these things, in this order:
%	* "mask_ch": numeric input that tells the function which channel has
%		the image data from which to derive the mask.  This can be used to
%		calculate the autocorrelation function for any arbitrarily-shaped
%		region, and perhaps even disjoint regions.  If passed as a scalar
%		numerical zero, or logical false, then the mask will not be used
%		and the autocorrelation function for the entire frame will be
%		calculated. Default, 2.
%		If this is not specified, but you still want to
%		specify other arguments, put empty brackets -- [] -- in place of
%		this argument.
%	* "yesplot": whether you want to plot the outcome.  Default, "false".  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	mask_ch = varargin{iArg}; else 
	mask_ch = 2;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	yesplot = varargin{iArg}; else
 	yesplot = false;
end%, iArg = iArg + 1;


vsep = strfind(filename,filesep);
vsep1 = vsep(end);
% vsep2 = vsep(end-1);
filenameshort = filename(vsep1+1:end);
pthname = filename(1:vsep1);

% Extract date using regular expression
expr = '\d{4}-\d{2}-\d{2}';  % pattern for YYYY-MM-DD
tokens = regexp(filename, expr, 'match');

% Assign to datum (use first match)
if ~isempty(tokens)
    datum = tokens{1};
else
    datum = '';
    warning('No date found in filename.');
end


[IM,lsminf1,lsminf2] = openczi(filename);
% lsminf1 = lsm.lsm;
% lsminf2 = lsminfo(filename);

scalings = lsminf2.VOXELSIZES*1e6; % microns per pixel.

% t = lsminf2.TimeStamps.TimeStamps;
% t = t - t(2);

r=bfGetReader(filename);
nSeries = r.getSeriesCount();
m = r.getMetadataStore();
nImages = m.getImageCount();
nPlanes = m.getPlaneCount(0);

t1 = zeros(nPlanes/2,1);
t2 = zeros(nPlanes/2,1);
for i=0:nPlanes-1
    
    if rem(i,2) == 1
        
        t1((i+1)/2,1) = m.getPlaneDeltaT(0,(i)).value;
         
    else
        
        t2((i+2)/2,1) = m.getPlaneDeltaT(0,(i)).value;
                    
    end
    
end

t = t1;
t = t-t(2);

% filenameshort = '2021-12-15_Cact-FRAP embryo 1 pos 1 z1';

% IM = lsm.data;
[H,W,num_ch,num_frames] = size(IM);


if num_frames < 5
	% skip this file if not enough frames
	data = 'Too few frames';
	return
end

%default
% data_ch = 2;
% mask_ch = 1;

% for 2022-10-11 etika frad czi data
data_ch = 1;
mask_ch = 2;

if num_ch >= 2
	I_dl = squeeze(IM(:,:,data_ch,:));
	I_nuc = squeeze(IM(:,:,mask_ch,:)); % red channel is best
else
	I_dl = squeeze(IM(:,:,data_ch,:));
	I_nuc = I_dl; % but some time series don't have it
end


if H < 50 || W < 50
	% skip this file if, after stabilizing, there is not enough image left
	data = 'Nuclei moved too much';
	return
end


% Call function to segment the nuclei.  The output of this function
% includes L, which is essentially a label array for every frame.

[L_nuc,L_cyt,nuc_data,cyt_data,Nucindex,side,mask_data,label] = ...
		find_nucleizoom_new_comp_set2(I_dl,I_nuc,t,scalings,datum,filenameshort,filenum,0);

% % Which nucleus is the one that is bleached? 
% [~,k_bleach] = max(nucdl(1,:)-nucdl(2,:));

% Enter the identity of the bleached nucleus
prompt = input('Please enter bleached nuc label: ');
idx_bleach = prompt; 
% idx_bleach = 16;


%=======================================================================
% use this to reconstruct and track a bleached nucleus if it's missing for
% some frames. For example used for filenum: 2, 14, 20

% uncomment to run the reconstruct the missing bleached nucleus

% if ismember(filenum, [2, 14, 20])
% 
%     % Enter the identity of the bleached nucleus
%     prompt = input('Please enter the number of frames to reconstruct the missing nucleus: ');
%     num_frames_fix = prompt; 
% 
%     [L_nuc,L_cyt,nuc_data,cyt_data,Nucindex] =...
%         reconstruct_missing_nucleus(mask_data,label,nuc_data,I_dl,I_nuc,t,...
%         idx_bleach,num_frames_fix,datum,filenameshort);
% 
%     % May need to reenter the identity of the bleached nucleus if nucleus was
%     % reconstructed
%     prompt = input('Please re-enter bleached nuc label: ');
%     idx_bleach = prompt;
% 
% end


%=======================================================================

% Extracting the nuclear and cytoplasmic concentrations as a
% function of time for each nucleus and cytoplasmic island.

[nucdl,cytdl] = dl_vs_t(I_dl,L_nuc,L_cyt,Nucindex);
n_nuc = size(nucdl,2);

%=======================================================================

% locs of surrounding nuclei -returns the locations of the closest neighbors. 
% Default number of nighbors set at 6 but can be changed

locs = persistent_closest_neighbors(nuc_data,idx_bleach,6); 

%=======================================================================
% plot and save the image results for the Neighbor nuc stats

[mean_oth_nuc,mean_oth_cyt] = plot_neighbors(nucdl,cytdl,idx_bleach,locs,t,datum,filenameshort,filenum);

%=======================================================================
% getting the nuclear and cytoplasmic intensities for bleached nuc/cyto...
% and the average nuc and cyto

y_nuc = nucdl(:,idx_bleach);
y_cyt = cytdl(:,idx_bleach);


% The rest of the nuclei and cytoplasm undergoes change over time.
%  How does that happen overall?

y_nucavg = meanDU(nucdl(:,setdiff(1:size(nucdl,2),idx_bleach)),2);
y_cytavg = meanDU(cytdl(:,setdiff(1:size(cytdl,2),idx_bleach)),2);


%
% Old variables:
%
L0 = NaN;
wsstats = NaN;

%
% metadata 2
%
data.filenameshort = filenameshort;
data.H = H;
data.W = W;
data.num_frames = num_frames;
data.num_ch = num_ch;
data.side = side;
data.nuc0 = L0;
data.wsstats = wsstats;
data.nuc_data = nuc_data;
data.cyt_data = cyt_data;
data.locs = locs;
data.t = t;
data.y_nuc = y_nuc;
data.y_cyt = y_cyt;
data.y_nucavg = y_nucavg;
data.y_cytavg = y_cytavg;
data.y_nucall = nucdl;
data.y_cytall = cytdl;
data.idx_bleach = idx_bleach;

%
% Metadata 3: placeholders
%
data.model = NaN;
data.yfit = NaN;
data.K_nuc = NaN;
data.k_out = NaN;
data.k_in = NaN;
data.c0 = NaN;
data.r2 = NaN;
data.delta1 = NaN;
data.delta2 = NaN;
data.covar = NaN;


%
% plotting invisibly, if asked for
%
if yesplot
	figure%('visib','off')
	set(gcf,'paperpositionmode','auto','pos',[20 355 1478 360]);
	
	subplot(1,3,1)
	imshow(I_dl(:,:,2),[])
	
	subplot(1,3,2)
	Ln = uint8(255*mDU(sum(L_nuc>0,3)));
	Lc = uint8(255*mDU(sum(L_cyt>0,3)));
	I1 = c1628(I_nuc(:,:,2));
	u = cat(3,I1,Lc,Ln);
	imshow(u)
	hold on
	
	
	subplot(1,3,3)
	plot(t(2:end),[y_nuc(2:end) y_cyt(2:end)],'o')
	xlabel('t_{all}')
	
	vsep = find(filename == filesep);
	vsep2 = vsep(end-1);
	filename1 = filename;
	filename1(vsep) = '/';
	
	%title(filename1(vsep2+1:end))
	
	saveas(gcf,[filename(1:end-4),'.fig'])
	close
end





