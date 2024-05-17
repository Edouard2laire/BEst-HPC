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
    folder_data = fullfile('/Users/edelaire1/Documents/Project/BEst-HPC/data', sProcess.options.token.Value  ,'out');
    if ~exist(folder_data)
        bst_error(' Unable to find the data');
        return; 
    end

    files = dir(fullfile(folder_data, '*.mat'));
    for iFile = 1:length(files)
        sData = load(fullfile(files(iFile).folder, files(iFile).name));
        bst_info = sData.bst_info;
        
        sStudy = bst_info.sStudy; 
        [sStudy, iStudy] = bst_get('Study', sStudy.FileName);

        
        %% Save results
        nData = length(sData.sResultsBst);
        bst_progress('text', 'Saving Results...');
        for iData=1:nData

            ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'results_MEM');
            ResultsMat = sData.sResultsBst(iData);
            ResultsMat = bst_history('add', ResultsMat, 'Import data', 'Import data from HPC');
            bst_save(ResultFile, ResultsMat, 'v7.3');

            db_add_data(iStudy, ResultFile, ResultsMat);
            OutputFiles{end+1} = ResultFile;
        end

    end
end


