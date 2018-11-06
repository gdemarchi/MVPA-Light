
%%% obob
addpath('/mnt/obob/staff/gdemarchi/git/obob_ownft/');
cfg = [];
cfg.package.svs = 'true';
cfg.package.obsolete = 'true';
obob_init_ft(cfg);

%%% the rest
addpath(genpath('/mnt/obob/staff/gdemarchi/git/MVPA-Light/'));
startup_MVPA_Light;

addpath /mnt/obob/staff/gdemarchi/mattools/;
addpath /mnt/obob/staff/gdemarchi/git/export_fig/;


rawDir= (['/mnt/obob/staff/gdemarchi/data/markov/raw/sss/']);
fileDir = (['/mnt/obob/staff/gdemarchi/data/markov/decoding/TG_trSNDteOM_prestim/final/']);

%%%%%%%%%%%%%%%%%%%%%%%  COMMON PART  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all;
clear tmpdata data % to stay on the safe side
conds={'random*','midminus*','midplus*','ordered*'};
trialinfos = [];

%% try another grad i.e. the non ica-ed
tmpFileGrad= dir([rawDir,'*',subJ,'_block01*']);
cur_fileGrad = [tmpFileGrad.folder,'/',tmpFileGrad.name];

%%% read the stuff
% not ica timegen_trRdSNDteALLOM_prestim_UNCLEANED_plusW_unbalanced_woSelfRepetitions_lda_Fs100_SelfRepNotRemoved
fileToRead  = dir([fileDir subJ '*UNCLEANED_plusW_ICAcleaned_unbalanced_woSelfRepetitions_lda_Fs100_reallyFinal*']);

curFile = fileToRead.name;
tmptg=load([fileDir curFile]);
TG_res  = tmptg.result_accTG_RdSNDPostStim_RdOM;
stuffForWeigths = tmptg.stuffForWeights_trRdSNDpost;
trainTime= tmptg.timeTrain;


%%% for each time point ...
for iTime =1:length(tmptg.timeTrain)
  tmppattern = mv_stat_activation_pattern(stuffForWeigths.cf{iTime}, stuffForWeigths.data(:,:,iTime), stuffForWeigths.trialinfo(:,1));
  weigths(:,iTime) =tmppattern(:,1); %channel x 3, here take the strongest SVD ... just to try ...
  fprintf('timepoint: %d done!\n', iTime);
end

%%% create the fake topography
clear tmpW
tmpW = [];
tmpW.avg = squeeze(weigths);
tmpW.time = trainTime; %fixin the saving ...
tmpW.grad= stuffForWeigths.grad;
tmpW.label= stuffForWeigths.label;
tmpW.dimord = 'chan_time';

% 
% %% check topos
figure;
cfg =[];
cfg.layout = 'neuromag306mag.lay';
ft_multiplotER(cfg,tmpW)

%%% quick source -init part

%% steal  headstuff ffrom marta
% /mnt/obob/staff/mpartyka/markov/SRC/data/head_aligned

headFile = ['/mnt/obob/staff/mpartyka/markov/SRC/data/head_aligned/*' subJ '*.mat'];
fName = dir(headFile);
load([fName.folder '/' fName.name]); %hdm mri_aligned mri_segmented shape

%%!!! make the grid point extenally selectable
%load mni_grid_1_5_cm_889pnts.mat 
%load mni_grid_1_cm_2982pnts.mat
load mni_grid_2_5_cm_191pnts.mat

%load standard_mri_better.mat
%% check meters all over the place !!!
cfg = [];
cfg.coordsys='neuromag';
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
%%%cfg.grid.resolution = 0.03;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri_aligned;
cfg.headmodel = hdm;
individual_grid               = ft_prepare_sourcemodel(cfg);

%%
stuffForWeigths.grad = ft_convert_units(stuffForWeigths.grad,'m');
%gradHipp = ft_convert_units(gradHipp,'m');

hdm  = ft_convert_units(hdm,'m');
individual_grid =  ft_convert_units(individual_grid,'m');

cfg=[];
cfg.channel = {'MEGMAG'};
cfg.vol=hdm;
cfg.grid=individual_grid;
%%cfg.grid.resolution = 0.03; %3 cm
cfg.grid.unit = 'm';
cfg.grad=stuffForWeigths.grad;
cfg.normalize='yes';
cfg.inwardshift = -0.2; %m?!
%cfg.moveinward = -2.5;
lf=ft_prepare_leadfield(cfg);
%cfg.grad=gradHipp;
%lfHipp=ft_prepare_leadfield(cfg);

%%  badly built input, RECHECK!!!!

nTrl = max(size(stuffForWeigths.data));
timeTrl = stuffForWeigths.time;

stuffForWeigths = rmfield(stuffForWeigths,'time');

for iTrl = 1:nTrl
  stuffForWeigths.trial{iTrl}   = squeeze(stuffForWeigths.data(iTrl,:,:));
  stuffForWeigths.time{iTrl}    = timeTrl;
end
%stuffForWeigths.trial ={stuffForWeigths.data};%keep an eye on this
stuffForWeigths.dimord = 'rpt_chan_time'; %maybe not neede in future
stuffForWeigths = rmfield(stuffForWeigths,'data');

%% do PCA to help the beaformer later 
%maybe here i shoul fix for grads as well
% or maybe not let's see ..
% cfg = [];
% cfg.method = 'pca';
% cfg.updatesens = 'no';
% cfg.channel = {'MEG*1'};
% comp = ft_componentanalysis(cfg, stuffForWeigths);
% 
% %%
% cfg = [];
% cfg.updatesens = 'no';
% cfg.component = comp.label(55:end); %through small comps away
% stuffForWeigths = ft_rejectcomponent(cfg, comp);

%%

cfg=[];
cfg.channel = {'MEGMAG'};
cfg.preproc.hpfilter   = 'yes';
cfg.preproc.hpfreq     = 1;
cfg.preproc.lpfilter   = 'yes';
cfg.preproc.lpfreq     = 45;
cfg.covariance         = 'yes';
%cfg.covariancewindow   = covariancewindow;

data_avg = ft_timelockanalysis(cfg, stuffForWeigths);

%% do lcmv beamforming (better: compute filter on averaged data)

% compute spatial filters
cfg=[];
cfg.method          = 'lcmv';
cfg.vol.unit        = 'm'; % th: dirty hack to make this work. as we always provide leadfields, we do not need a vol...
cfg.grid            = lf;
cfg.projectnoise    = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '15%'; % 5% doesnt' seem to be enough according to Litvak then 15%
cfg.lcmv.powmethod = 'yes';
lcmvall=ft_sourceanalysis(cfg, data_avg);
%cfg.grid            = lfHipp;
%lcmvallHipp=ft_sourceanalysis(cfg, data_avg);

% normal
beamfilts = cat(1,lcmvall.avg.filter{:});

data_source = stuffForWeigths;
data_source=rmfield(data_source,'trial');

%%% i need to project the weights only, i.e. the topography
data_source.avg = beamfilts*tmpW.avg;

data_source.label = cellstr(num2str([1:sum(lf.inside)]'));
data_source.dimord = 'chan_time'; %maybe not neede in future
data_source.time = tmpW.time;

% hipp
% beamfiltsHipp = cat(1,lcmvallHipp.avg.filter{:});
% 
% data_sourceHipp = stuffForWeigths;
% data_sourceHipp=rmfield(data_sourceHipp,'trial');
% data_sourceHipp.avg = beamfiltsHipp*tmpW.avg;
% data_sourceHipp.grad = gradHipp;
% data_sourceHipp.label = cellstr(num2str([1:sum(lfHipp.inside)]'));
% data_sourceHipp.dimord = 'chan_time'; %maybe not neede in future
% data_sourceHipp.time = tmpW.time;

%% here I could save, since the rest ist swiftly doable later ..

outFile = [ subJ '_weights_beamed_15pcRegFac_yesICA_25cmgrid_reallyFinal.mat'];
save (fullfile(outDir, outFile),'data_source*','tmpW*','-v7.3');

% %%
% figure;
% plot(data_source.time, abs(mean(data_source.avg)))
% %plot(data_source.time, abs(mean(data_source.avg(individual_grid.inside))))
% 
% %%
% cfg=[];
% cfg.baseline=[.5 .55];
% cfg.baselinetype='absolute';
% data_sourcebl=obob_svs_timelockbaseline(cfg, data_source);
% 
% %% vs2s
% 
% cfg=[];
% cfg.sourcegrid = template_grid%%individual_grid% 
% cfg.parameter={'avg'};
% cfg.toilim=[.09 .15];
% cfg.mri = mri% mri_aligned; %template mri %load standard_mri_better.mat
% source2plot = obob_svs_virtualsens2source(cfg, data_source);
% %source2plot.pow = abs(source2plot.avg).*source2plot.inside;
% 
% 
% %% single
% figure;
% cfg = [];
% cfg.funparameter = 'pow';
% %cfg.maskparameter = cfg.funparameter %'mask2';
% %cfg.atlas = '/Users/gianpaolo/git/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii';
% %cfg.roi = {'Calcarine_L', 'Calcarine_R'}
% %cfg.interactive = 'no';
% %cfg.crosshair = 'no';%'yes';
% %cfg.axis = 'off';
% %cfg.colorbar = 'yes';
% cfg.funcolormap = 'hot';
% cfg.funcolorlim   = 'maxabs'%[0 1].*10^27%'zeromax' ;%'maxabs'%[-1 1].*10^22;
% %cfg.funcolorlim   = [0 5];%
% %cfg.location      = [173 46 138];
% %cfg.locationcoordinates = 'voxel';
% ft_sourceplot(cfg, source2plot);
% 
% %at the single subject level they seem to be fine ...
% %let's take a look tomorrow on all the subjects
% 
% % no difference between with and without Hipp and PCA
% % no difference on the GA average with and without ICA
% 
% % the only difference seems to be in the topographies i.e. weigths after
% % ica
% 
% 

