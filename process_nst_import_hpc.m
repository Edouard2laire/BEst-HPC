function varargout = process_nst_import_hpc( varargin )
% process_nst_import_hpc:  import all the results comming from the
% HPC-computers

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
    sProcess.Comment     = 'Import data from HPC';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Simulation'};
    sProcess.Index       = 3005;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.token.Comment = 'Token';
    sProcess.options.token.Type    = 'text';
    sProcess.options.token.Value   =  '';
end
%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    OutputFiles = {};
    folder_data = fullfile('/Users/edelaire1/Documents/Project/wMEM-fnirs', sProcess.options.token.Value  ,'out');
    if ~exist(folder_data)
        bst_error(' Unable to find the data');
        return; 
    end

    files = dir(fullfile(folder_data, '*.mat'));
    for iFile = 1:length(files)
        sData = load(fullfile(files(iFile).folder, files(iFile).name));
        
        sStudy = sData.sStudy; 
        [sStudy, iStudy] = bst_get('Study', sStudy.FileName);

        ChannelMat = in_bst_channel(sStudy.Channel.FileName);
        headmodel = in_bst_headmodel(sStudy.HeadModel.FileName);
        headmodel.FileName = sStudy.HeadModel.FileName;

        sInputs = struct('FileName', sData.OPTIONS.DataFile,'iStudy', iStudy);

        if ~isfield(sData, 'Hb_sources')
            bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');
            % Compute dHb
            hb_extinctions = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths);
            hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L
        
            Hb_sources = zeros(size(sData.dOD_sources,1), 3, size(sData.dOD_sources,3));
            valid_nodes = sData.diagnosis.valid_nodes;
            for idx=1:length(valid_nodes)
                inode = valid_nodes(idx);
                Hb_sources(inode, 1:2, :) = pinv(hb_extinctions) * ...
                                            squeeze(sData.dOD_sources(inode, :, :));
            
            end
            Hb_sources(:,3,:) = squeeze(sum(Hb_sources, 2));
            sData.Hb_sources = Hb_sources;
        end

        %% Save results
        bst_progress('text', 'Saving Results...');
        for iwl=1:1:length(ChannelMat.Nirs.Wavelengths)
             swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];

             [sStudy, ResultFile] = process_nst_cmem('add_surf_data',squeeze(sData.dOD_sources(:,iwl,:)), sData.OPTIONS.DataTime, headmodel, ...
                                                  [sData.diagnosis(iwl).Comment     ' | ' swl 'nm'], ...
                                                   sInputs, sStudy, sData.diagnosis(iwl).Comment, ...
                                                   'OD', 0,1, sData.diagnosis(iwl));
             OutputFiles{end+1} = ResultFile;
        end
        
        hb_unit_factor = 1e6;
        hb_unit = '\mumol.l-1';
        hb_types = {'HbO', 'HbR','HbT'};
        for ihb=1:3
            [sStudy, ResultFile] = process_nst_cmem('add_surf_data',squeeze(sData.Hb_sources(:,ihb,:)) .* hb_unit_factor,...
                                                 sData.OPTIONS.DataTime, headmodel, ...
                                                 [sData.diagnosis(iwl).Comment     ' | ' hb_types{ihb}], ...
                                                 sInputs, sStudy, sData.diagnosis(iwl).Comment, ...
                                                 hb_unit, 0);    
            OutputFiles{end+1} = ResultFile;
        end

    end
end


