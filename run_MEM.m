function run_MEM(data, options_file)
    rng("shuffle")
    
    fprintf('Using ')
    system('hostname');
    disp('')
    fprintf('Data file : %s \n', data)
    fprintf('Option file : %s \n', options_file)


    addpath(genpath('best-brainstorm-master/best'));
        
    input_struct = load(data);
    sData        = input_struct.sData;
    bst_info     = input_struct.bst_info;

    fprintf('Data  : %s \n', bst_info.Comment)
    fprintf('--------------------------------\n')

    output_file     = strrep(data,'in','out');
    output_folder   = fileparts(output_file);
    if ~exist(output_folder)
        mkdir(output_folder)
    end
    
    % Loading options 
    OPTIONS = jsondecode(fileread(options_file));
    OPTIONS = OPTIONS.MEMpaneloptions;
    OPTIONS.automatic.stand_alone =   1;
    OPTIONS.automatic.Comment     = 'MEM';
    OPTIONS.automatic.selected_samples = [];

    if isfield(OPTIONS,'wavelet') && isfield(OPTIONS.wavelet,'selected_scales') && size(OPTIONS.wavelet.selected_scales,1) > 1
        OPTIONS.wavelet.selected_scales = OPTIONS.wavelet.selected_scales';
    end

    nb_nodes = size(bst_info.sCortex.Vertices, 1);
    nb_samples = length(sData(1).OPTIONS.mandatory.DataTime);
    nb_data  = length(sData);

    
    % try to start parralel port
    % if isempty(gcp('nocreate'))
    %     local_cluster = parcluster('local');
    %     parpool(local_cluster, 12);
    % end
    
    sources = zeros(nb_nodes, nb_data, nb_samples);
    sResultsBst = repmat(getTemplate(),1,nb_data);
    for iData=1:nb_data
        fprintf('MEM > Computing MEM for data %d\n', iData); 
        
        sDataWl     = sData(iData);
        MEMoptions  = be_struct_copy_fields( sDataWl.OPTIONS, OPTIONS, [],1 );

        %% launch MEM 
        [Results, O_updated] = be_main_call(sDataWl.HM, MEMoptions);

        % MEM results
        grid_amp = zeros(nb_nodes, nb_samples); 
        grid_amp(bst_info.valid_nodes,:) = Results.ImageGridAmp;
        
        sResultsBst(iData).ImageGridAmp = grid_amp;
        sResultsBst(iData).Time = sDataWl.OPTIONS.mandatory.DataTime;
        sResultsBst(iData).DataFile = bst_info.DataFile;
        sResultsBst(iData).HeadModelFile = bst_info.HeadModelFile;
        sResultsBst(iData).SurfaceFile = bst_info.SurfaceFile;
        
        if isfield(bst_info,'hb_extinctions')
            sResultsBst(iData).Comment      = [O_updated.Comment     ' | '  'WL' num2str( bst_info.Wavelengths(iData)) 'nm'];
        else
            sResultsBst(iData).Comment      = O_updated.Comment;
        end
        sResultsBst(iData).Options      = O_updated;
        if isfield(bst_info,'hb_extinctions')
             sResultsBst(iData).DisplayUnits = 'OD';
        end
        sources(:, iData, :)  = grid_amp;
    end

    if isfield(bst_info,'hb_extinctions')
        hb_extinctions = bst_info.hb_extinctions;
    
        Hb_sources = zeros(nb_nodes, 3, nb_samples);
        for idx=1:length(bst_info.valid_nodes)
            inode = bst_info.valid_nodes(idx);
            Hb_sources(inode, 1:2, :) = pinv(hb_extinctions) * ...
                                        squeeze(sources(inode, :, :));
        end
        Hb_sources(:,3,:) = squeeze(sum(Hb_sources, 2));

        sResultsBstHb = repmat(getTemplate(),1,3);
        hb_unit_factor = 1e6;
        hb_unit = '\mumol.l-1';
        hb_types = {'HbO', 'HbR','HbT'};

        for ihb=1:3
            sResultsBstHb(ihb).ImageGridAmp = squeeze(Hb_sources(:,ihb,:)) .* hb_unit_factor;
            sResultsBstHb(ihb).Time = sDataWl.OPTIONS.mandatory.DataTime;
            sResultsBstHb(ihb).DataFile = bst_info.DataFile;
            sResultsBstHb(ihb).HeadModelFile = bst_info.HeadModelFile;
            sResultsBstHb(ihb).SurfaceFile = bst_info.SurfaceFile;
            sResultsBstHb(ihb).Comment      = [O_updated.Comment    ' | ' hb_types{ihb}] ;
            sResultsBstHb(ihb).Options      = O_updated;
            sResultsBstHb(ihb).DisplayUnits = hb_unit;
        end
        sResultsBst = [sResultsBst , sResultsBstHb];

    end

    %delete(gcp('nocreate'))
    save(output_file, 'bst_info', 'OPTIONS', 'sResultsBst','-v7.3');

end


function template = getTemplate()
        template = struct('ImagingKernel', [], ...
                          'ImageGridAmp',  [], ...
                          'Std',           [], ...
                          'Whitener',      [], ...
                          'SourceDecompSa',[], ...
                          'SourceDecompVa',[], ...
                          'nComponents',   1, ...
                          'Comment',       '', ...
                          'Function',      '', ...
                          'Time',          [], ...
                          'DataFile',      '', ...
                          'HeadModelFile', '', ...
                          'HeadModelType', 'surface', ...
                          'ChannelFlag',   [], ...
                          'GoodChannel',   [], ...
                          'SurfaceFile',   [], ...
                          'Atlas',         [], ...
                          'GridLoc',       [], ...
                          'GridOrient',    [], ...
                          'GridAtlas',     [], ...
                          'Options',       [], ...
                          'ColormapType',  [], ...
                          'DisplayUnits',  [], ...
                          'ZScore',        [], ...
                          'nAvg',          [], ...
                          'Leff',          [], ...
                          'History',       []);
end

