function varargout = process_nst_export_hpc( varargin )
% process_nst_export_hpc:  Compute source-localization using Bootstrap
% and MEM

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
% Authors: Edouard Delaire, 2022

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

    sProcess.options.auto_neighborhood_order.Comment = 'Set neighborhood order automatically (default)';
    sProcess.options.auto_neighborhood_order.Type    = 'checkbox';
    sProcess.options.auto_neighborhood_order.Value   = 1;

    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {3, 'cm',2};

end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    OutputFiles = {};

    token = char(floor(26*rand(1, 10)) + 65); 
    script_path = '/Users/edelaire1/Documents/Project/wMEM-fnirs';

    folder_out = fullfile('/Users/edelaire1/Documents/Project/wMEM-fnirs/data', token );
    mkdir(fullfile(folder_out,'in'));

    
    %% Load head model
    sStudy = bst_get('Study', sInputs.iStudy);
    if isempty(sStudy.iHeadModel)
        bst_error('No head model found. Consider running "NIRS -> Compute head model"');
        return;
    end
    
    nirs_head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
    nirs_head_model.FileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    
    cortex = in_tess_bst(nirs_head_model.SurfaceFile);
    
    fID = fopen(fullfile(script_path, sprintf('launch_script_%s.sh', token)), 'w+');

    for iInput = 1:length(sInputs)
        sData = struct('HM', struct(), 'OPTIONS', struct());


        if strcmp(sInputs(iInput).FileType, 'data')     % Imported data structure
            sDataIn = in_bst_data(sInputs(iInput).FileName);
        elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
            sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
        end
    
        ChannelMat = in_bst_channel(sInputs(iInput).ChannelFile);
        if ~isfield(ChannelMat.Nirs, 'Wavelengths')
            bst_error(['cMEM source reconstruction works only for dOD data ' ... 
                       ' (eg do not use MBLL prior to this process)']);
            return;
        end

        thresh_dis2cortex       =  sProcess.options.thresh_dis2cortex.Value{1}  / 100;
        valid_nodes             = nst_headmodel_get_FOV(ChannelMat, cortex, thresh_dis2cortex,sDataIn.ChannelFlag );

        fprintf('MEM > Estimating neighborhood order\n'); 
        if sProcess.options.auto_neighborhood_order.Value
            swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(1))];
            n_channel = sum(strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)');
        
            nbo = process_nst_cmem('estimate_nbo',cortex, valid_nodes, n_channel, 1 );
            fprintf('MEM > Using a NBO of %d\n', nbo); 
        end

        nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);
        for iwl=1:nb_wavelengths
            swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
            selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
            

            % Define Head model
            gain = nst_headmodel_get_gains(nirs_head_model,iwl, ChannelMat.Channel, find(selected_chans));
        
            % Remove 0 from the gain matrixHeadModel.Gain(1).matrix
            tmp = gain(:,valid_nodes);
            tmp(tmp == 0) = min(tmp(tmp > 0));
            
            HM = struct();
            HM.vertex_connectivity = cortex.VertConn(valid_nodes, valid_nodes);
            HM.Gain(1).matrix = tmp;
            HM.Gain(1).modality = 'NIRS';
            sData(iwl).HM = HM;


            % Define data 
            OPTIONS = struct(); 
            OPTIONS.mandatory.DataTime          = round(sDataIn.Time,6);
            OPTIONS.mandatory.Data              = sDataIn.F(selected_chans,:);
            OPTIONS.mandatory.DataTypes         =  {ChannelMat.Channel(1).Type};
            OPTIONS.mandatory.ChannelTypes      = {ChannelMat.Channel(selected_chans).Type};
            OPTIONS.mandatory.GoodChannel       = ones(sum(selected_chans), 1);
            OPTIONS.mandatory.ChannelFlag       = ones(sum(selected_chans), 1);
            
            sData(iwl).OPTIONS = OPTIONS;
        end



        bst_info                 = struct();
        bst_info.Comment         = sInputs(iInput).Comment;
        bst_info.sStudy          = sStudy;
        bst_info.DataFile        = sInputs(iInput).FileName;
        bst_info.SurfaceFile     = nirs_head_model.SurfaceFile;
        bst_info.HeadModelFile   = nirs_head_model.FileName;
        bst_info.sCortex         = cortex;
        bst_info.nbo             = nbo; 
        bst_info.valid_nodes     = valid_nodes;
        bst_info.hb_extinctions  = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths)./10;% mm-1.mole-1.L

        sOutput = struct(  'sData', sData, ...
                           'bst_info',bst_info); 
    
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


