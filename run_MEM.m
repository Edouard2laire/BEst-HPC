function run_MEM(data, options_file)
    rng("shuffle")
    
    fprintf('Using ')
    system('hostname');
    disp('')
    fprintf('Data file : %s \n', data)
    fprintf('Option file : %s \n', options_file)


    addpath(genpath('best-brainstorm-master/best'));
        
    input_struct = load(data);
    sData = input_struct.sData;
    bst_info = input_struct.bst_info;

    fprintf('Data  : %s \n', bst_info.Comment)
    fprintf('--------------------------------\n')

    output_file = strrep(data,'in','out');
    output_folder = fileparts(output_file);
    if ~exist(output_folder)
        mkdir(output_folder)
    end
    
    % Loading options 
    OPTIONS = jsondecode(fileread(options_file));
    OPTIONS = OPTIONS.MEMpaneloptions;
    OPTIONS.automatic.stand_alone =   1;
    OPTIONS.automatic.Comment     = 'MEM';
    OPTIONS.wavelet.selected_scales = OPTIONS.wavelet.selected_scales';

    nb_nodes = size(bst_info.sCortex.Vertices, 1);
    nb_samples = length(sData(1).OPTIONS.mandatory.DataTime);
    nb_wavelengths  = length(sData);

    
    % try to start parralel port
    % if isempty(gcp('nocreate'))
    %     local_cluster = parcluster('Processes');
    %     parpool(local_cluster);
    % end

    dOD_sources = zeros(nb_nodes, nb_wavelengths, nb_samples);
    diagnosis   = [];

    for iwl=1:nb_wavelengths
        fprintf('MEM > Computing MEM for wavelength %d\n', iwl); 
        
        sDataWl = sData(iwl);
        

        MEMoptions = be_struct_copy_fields( sDataWl.OPTIONS, OPTIONS, [],1 );

        %% launch MEM (cMEM only in current version)
        [Results, O_updated] = be_main_call(sDataWl.HM, MEMoptions);

        %cMEM results
        grid_amp = zeros(nb_nodes, nb_samples); 
        grid_amp(bst_info.valid_nodes,:) = Results.ImageGridAmp;
        
        dOD_sources(:, iwl, :)  = grid_amp;
        diagnosis          = [diagnosis Results.MEMoptions.automatic];
    end

    hb_extinctions = bst_info.hb_extinctions;

    Hb_sources = zeros(nb_nodes, 3, nb_samples);
    for idx=1:length(bst_info.valid_nodes)
        inode = bst_info.valid_nodes(idx);
        Hb_sources(inode, 1:2, :) = pinv(hb_extinctions) * ...
                                    squeeze(dOD_sources(inode, :, :));
    
    end
    Hb_sources(:,3,:) = squeeze(sum(Hb_sources, 2));

    %delete(gcp('nocreate'))
    %save(output_file, 'bst_info', 'OPTIONS', 'dOD_sources','Hb_sources', 'diagnosis','-v7.3');
end

