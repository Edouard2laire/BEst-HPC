function varargout = process_nst_export_hpc( varargin )
% process_nst_export_hpc: Export data to compute MEM using the standalone
% version on HPC computers

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2018 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Edouard Delaire, 2022 - 2024

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Export data to HPC';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Simulation'};
    sProcess.Index       = 3004;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.auto_neighborhood_order.Comment = '(NIRS only) Set neighborhood order automatically (default)';
    sProcess.options.auto_neighborhood_order.Type    = 'checkbox';
    sProcess.options.auto_neighborhood_order.Value   = 1;

    sProcess.options.thresh_dis2cortex.Comment = '(NIRS only)  Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {3, 'cm',2};

        % Option: Sensors selection
    sProcess.options.sensortypes.Comment = 'Sensor types:&nbsp;&nbsp;&nbsp;&nbsp;';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG, NIRS';

end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)

    OutputFiles = {};

    token = char(floor(26*rand(1, 10)) + 65); 
    script_path = '/Users/edelaire1/Documents/Project/BEst-HPC/data';

    folder_out = fullfile(script_path, token );
    mkdir(fullfile(folder_out,'in'));

    
    %% Load head model
    sStudy = bst_get('Study', sInputs.iStudy);
    if isempty(sStudy.iHeadModel)
        bst_error('No head model found.');
        return;
    end

    HeadModel = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
    HeadModel.FileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    
    sCortex = in_tess_bst(HeadModel.SurfaceFile);
    
    fID = fopen(fullfile(script_path, sprintf('launch_script_%s.sh', token)), 'w+');

    for iInput = 1:length(sInputs)
        sData = struct('HM', struct(), 'OPTIONS', struct());

        if strcmp(sInputs(iInput).FileType, 'data')     % Imported data structure
            sDataIn = in_bst_data(sInputs(iInput).FileName);
        elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
            sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
        end
    
        ChannelMat = in_bst_channel(sInputs(iInput).ChannelFile);
        isNIRS = isfield(ChannelMat,'Nirs');

        if isNIRS
            thresh_dis2cortex       =  sProcess.options.thresh_dis2cortex.Value{1}  / 100;
            valid_nodes             =  nst_headmodel_get_FOV(ChannelMat, sCortex, thresh_dis2cortex,sDataIn.ChannelFlag );
            nb_wavelengths          = length(ChannelMat.Nirs.Wavelengths);
        
            fprintf('MEM > Estimating neighborhood order\n'); 
                
            swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(1))];
            n_channel = sum(strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)');
            
            nbo = process_nst_cmem('estimate_nbo',sCortex, valid_nodes, n_channel, 1 );
            fprintf('MEM > Using a NBO of %d\n', nbo); 
        else
            valid_nodes             = 1:size(sCortex.Vertices,1);
            nb_wavelengths          = 1;
            nbo                     = 4;
        end


        for iwl=1:nb_wavelengths

            if isNIRS
                swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
                GoodChannel = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
            
                % Define Head model
                tmp = nst_headmodel_get_gains(HeadModel,iwl, ChannelMat.Channel, find(GoodChannel));
            
                % Remove 0 from the gain matrixHeadModel.Gain(1).matrix
                Gain = tmp(:,valid_nodes);
                Gain(Gain == 0) = min(Gain(Gain > 0));
            else
                GoodChannel = good_channel(ChannelMat.Channel, sDataIn.ChannelFlag, sProcess.options.sensortypes.Value);
                % Apply current SSP projectors
                if ~isempty(ChannelMat.Projector)
                    % Rebuild projector in the expanded form (I-UUt)
                    Proj = process_ssp2('BuildProjector', ChannelMat.Projector, [1 2]);
                    % Apply projectors
                    if ~isempty(Proj)
                        % Get all sensors for which the gain matrix was successfully computed
                        iGainSensors = find(sum(isnan(HeadModel.Gain), 2) == 0);
                        % Apply projectors to gain matrix
                        HeadModel.Gain(iGainSensors,:) = Proj(iGainSensors,iGainSensors) * HeadModel.Gain(iGainSensors,:);
                    end
                end
                % Select only good channels
                Gain = HeadModel.Gain(GoodChannel, :);

                if any(ismember(unique({ChannelMat.Channel(GoodChannel).Type}), {'EEG','ECOG','SEEG'}))
                    % Create average reference montage
                    sMontage = panel_montage('GetMontageAvgRef', [], ChannelMat.Channel(GoodChannel), ChannelFlag(GoodChannel), 0);
                    Gain = sMontage.Matrix * Gain;
               end

                Gain = bst_gain_orient(Gain, HeadModel.GridOrient);
            end

            HM = struct();
            HM.vertex_connectivity = sCortex.VertConn(valid_nodes, valid_nodes);
            HM.Gain(1).matrix = Gain;
            HM.Gain(1).modality = ChannelMat.Channel(GoodChannel(1)).Type;
            sData(iwl).HM = HM;


            % Define data 
            OPTIONS = struct(); 
            OPTIONS.mandatory.DataTime          = round(sDataIn.Time,6);
            OPTIONS.mandatory.Data              = sDataIn.F(GoodChannel,:);
            OPTIONS.mandatory.DataTypes         =  {ChannelMat.Channel(GoodChannel(1)).Type};
            OPTIONS.mandatory.ChannelTypes      = {ChannelMat.Channel(GoodChannel).Type};
            OPTIONS.mandatory.GoodChannel       = GoodChannel;
            OPTIONS.mandatory.ChannelFlag       = ones(length(GoodChannel), 1);
            
            sData(iwl).OPTIONS = OPTIONS;
        end

        bst_info                 = struct();
        bst_info.Comment         = sInputs(iInput).Comment;
        bst_info.sStudy          = sStudy;
        bst_info.DataFile        = sInputs(iInput).FileName;
        bst_info.SurfaceFile     = HeadModel.SurfaceFile;
        bst_info.HeadModelFile   = HeadModel.FileName;
        bst_info.sCortex         = sCortex;
        bst_info.nbo             = nbo; 
        bst_info.valid_nodes     = valid_nodes;
        if isNIRS
            bst_info.hb_extinctions  = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths)./10;% mm-1.mole-1.L
        end

        %% Run MEM
        [~,sOutput_name] = fileparts(sInputs(iInput).FileName);
        save(fullfile(folder_out, 'in', sprintf('%s.mat', sOutput_name)), 'sData','bst_info');

        fprintf(fID, 'qsub -j y -o logs/%s.txt -pe smp 16 -S /bin/bash  -cwd -q matlab.q -N MEM_%s ./start_hpc.sh %s %s \n', ...
                       sprintf('%d_%s',iInput,token), sprintf('%d_%s',iInput,token), fullfile(token, 'in', sprintf('%s.mat', sOutput_name)),  'wMEM_options.json');
    end
    fclose(fID);

     fID = fopen(fullfile(folder_out, 'README.txt'), 'w+');
     fprintf(fID, 'Simulation created on %s \n', char(datetime('today')));
     fprintf(fID, 'To tranfert the data to concordia, execute\n');
     fprintf(fID, 'rsync --times --progress --update --recursive ~/Documents/Project/wMEM-fnirs/%s mcgill:/mfip/mfip1/edelaire/HPC/wMEM-fnirs  \n', token);
     fprintf(fID, 'rsync --times --progress --update --recursive ~/Documents/Project/wMEM-fnirs/%s mcgill:/mfip/mfip1/edelaire/HPC/wMEM-fnirs  \n', sprintf('launch_script_%s.sh', token));

    fprintf(fID, 'To launch the script, execute\n');
    fprintf(fID, 'cd /mfip/mfip1/edelaire/HPC/wMEM-fnirs \n');
    fprintf(fID, '%s\n', fullfile(script_path, sprintf('launch_script_%s.sh', token)));

    fprintf(fID, 'To collect the data from concordia, execute\n');
    fprintf(fID, 'rsync -pE --times --progress --update --recursive edelaire@perf-imglab07:/NAS/home/edelaire/Documents/Project/wMEM-fnirs/%s  ~/Documents/Project/wMEM-fnirs \n', token);
    fprintf(fID, 'rsync -pE --times --progress --update --recursive edelaire@perf-imglab07:/NAS/home/edelaire/Documents/Project/wMEM-fnirs/logs  ~/Documents/Project/wMEM-fnirs \n');

    fclose(fID);

end


