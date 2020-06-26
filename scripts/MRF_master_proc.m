%% process MRF data

%% dependencies

% should already be added in batch_proc.m

%% params

%% load k-space trajectories

load([dir_in fn_ksp])

%% load MRF raw data

load([dir_in fn_MRF_raw])

if exist('skiprecon','var') && skiprecon == 0
    %% reconstruct and coil combine raw data

    input_img_recon.img_data = output_MRF_raw.data;
    input_img_recon.output_ksp_traj = output_ksp_traj;
    input_img_recon.ecalib_threshold = 0.0; %0.001; % sensitivity threshold

    output_img_recon = MRF_img_recon_bart( input_img_recon );

    %% load B1 map

    load([dir_in fn_B1])
    output_img_recon.B1_map = output_B1_data.B1_map;

    % output_img_recon.B1_map = ones(input_img_recon.effMtx); % test

    %% save image recon

    save([dir_out fn_MRF_img_recon],'output_img_recon')
    
else
    load([dir_out fn_MRF_img_recon])
end


%% make or load dictionary

TE_s = output_MRF_raw.params.TE_s;% s
nomFlip_deg = output_MRF_raw.params.nomFlip_deg; % nominal flip angle in degrees
TRbase_s = output_MRF_raw.params.TRbase_s; % s

if doMagdict_epg == 1

    % load acquisition vectors
    data_struct = importdata([dir_in fn_MRF_seq_params]);

    % set dictionary construction parameters
    input_dict.FA_v = data_struct.data(:,1)*nomFlip_deg; % deg
    input_dict.phi_v = data_struct.data(:,2); % deg
    input_dict.TR_v = TRbase_s*1000 + data_struct.data(:,3); % ms
    input_dict.TE_v = TE_s*1000 + data_struct.data(:,4); % ms
    if ndims( output_MRF_raw.data ) == 2
        input_dict.nreps = size( output_MRF_raw.data, 1 );
    elseif ndims(  output_MRF_raw.data ) == 3
        input_dict.nreps = size( output_MRF_raw.data, 2 );
    end

    input_dict.delk = 1; % step between states equal to a full dephasing imparted by crusher gradient
    input_dict.szomega = 101; % number of factors of k to include in phase history
    input_dict.reduce = 1;

    % plot sequence parameters
    figure(1); clf;
    subplot(411)
    plot(input_dict.FA_v); ylabel('degrees'); title(['FA, TI is ' num2str(input_dict.TI) ' ms'])
    subplot(412)
    plot(input_dict.phi_v); ylabel('degrees'); title('phase')
    subplot(413)
    plot(input_dict.TR_v); ylabel('msec'); title('TR')
    subplot(414)
    plot(input_dict.TE_v); ylabel('msec'); title('TE')
    drawnow

    % do dictionary construction        
    output_dict_epg = MRF_dict_B1(input_dict);
    
    % save result
    output_dict_epg.fn = fn_MRF_water_dict_epg;
    save([dir_out fn_MRF_water_dict_epg],'output_dict_epg')
    
else
 
    load([dir_out fn_MRF_water_dict_epg]);
    
end

if doMagdict_pepg == 1
    
     % load acquisition vectors
    data_struct = importdata([dir_in fn_MRF_seq_params]);

    % set dictionary construction parameters
    input_dict.FA_v = data_struct.data(:,1)*nomFlip_deg; % deg
    input_dict.phi_v = data_struct.data(:,2); % deg
    input_dict.TR_v = TRbase_s*1000 + data_struct.data(:,3); % ms
    input_dict.TE_v = TE_s*1000 + data_struct.data(:,4); % ms
    if ndims( output_MRF_raw.data ) == 2
        input_dict.nreps = size( output_MRF_raw.data, 1 );
    elseif ndims(  output_MRF_raw.data ) == 3
        input_dict.nreps = size( output_MRF_raw.data, 2 );
    end
    
    % get rf_shape for ssEPG
    rf_shape = get_rf_shape( shape_in.fn,shape_in.N_t,nomFlip_deg,shape_in.sl_thick_m,shape_in.tau_s,shape_in.n_cycles_per_crush,shape_in.N_res,shape_in.fov_factor,1 );
    rf_shape.shape_in = shape_in;
    
    input_dict.np = 50;
    input_dict.delk = 1; % step between states equal to a full dephasing imparted by crusher gradient
    input_dict.szomega = 101; % number of factors of k to include in phase history
    input_dict.rf_shape = rf_shape;
    input_dict.reduce = 1;
    
    output_dict_pepg = MRF_dict_generator_pEPG(input_dict);    
    
    % save result
    output_dict_pepg.fn = fn_MRF_water_dict_pepg;
    save([dir_out fn_MRF_water_dict_ssepg],'output_dict_ssepg')
    
else
    
    load([dir_out fn_MRF_water_dict_pepg]);
    output_dict_pepg = output_dict;
    output_dict_pepg.fn = fn_MRF_water_dict_pepg;
    clear output_dict;
    
end
    
if doMagdict_ssepg == 1
    
    % load acquisition vectors
    data_struct = importdata([dir_in fn_MRF_seq_params]);

    % set dictionary construction parameters
    input_dict.FA_v = data_struct.data(:,1)*nomFlip_deg; % deg
    input_dict.phi_v = data_struct.data(:,2); % deg
    input_dict.TR_v = TRbase_s*1000 + data_struct.data(:,3); % ms
    input_dict.TE_v = TE_s*1000 + data_struct.data(:,4); % ms
    if ndims( output_MRF_raw.data ) == 2
        input_dict.nreps = size( output_MRF_raw.data, 1 );
    elseif ndims(  output_MRF_raw.data ) == 3
        input_dict.nreps = size( output_MRF_raw.data, 2 );
    end
    
    % get rf_shape for ssEPG
    rf_shape = get_rf_shape( shape_in.fn,shape_in.N_t,nomFlip_deg,shape_in.sl_thick_m,shape_in.tau_s,shape_in.n_cycles_per_crush,shape_in.N_res,shape_in.fov_factor,1 );
    
    input_dict.delk = rf_shape.N_cycle;
    input_dict.szomega =  rf_shape.Q;
    input_dict.rf_shape = rf_shape;
    input_dict.reduce = 1;
    
    output_dict_ssepg = MRF_dict_generator_ssEPG_par(input_dict);    
    
    % save result
    output_dict_ssepg.fn = fn_MRF_water_dict_ssepg;
    save([dir_out fn_MRF_water_dict_ssepg],'output_dict_ssepg')

else

    load([dir_out fn_MRF_water_dict_ssepg]);
    output_dict_ssepg = output_dict;
    output_dict_ssepg.fn = fn_MRF_water_dict_ssepg;
    clear output_dict;

end

%% do a standard dictionary match if applicable

if doTmaps == 1

    reduce_flag = 1;

    output_MRF_match_epg = MRF_dict_match_B1( output_img_recon, output_dict_epg, reduce_flag );

    figure(1000)
    imagesc(output_MRF_match_epg.T1_map)
    axis image
    colorbar()
    title('T1 map')

    figure(1001)
    imagesc(output_MRF_match_epg.T2_map)
    axis image
    colorbar()
    title('T2 map')

    figure(1002)
    imagesc(abs(output_MRF_match_epg.M0_map))
    axis image
    colorbar()
    title('Magnitude of M0')
    
    output_MRF_match_pepg = MRF_dict_match_B1( output_img_recon, output_dict_pepg, reduce_flag );

    figure(2000)
    imagesc(output_MRF_match_pepg.T1_map)
    axis image
    colorbar()
    title('T1 map')

    figure(2001)
    imagesc(output_MRF_match_pepg.T2_map)
    axis image
    colorbar()
    title('T2 map')

    figure(2002)
    imagesc(abs(output_MRF_match_pepg.M0_map))
    axis image
    colorbar()
    title('Magnitude of M0')
    
    output_MRF_match_ssepg = MRF_dict_match_B1( output_img_recon, output_dict_ssepg, reduce_flag );

    figure(3000)
    imagesc(output_MRF_match_ssepg.T1_map)
    axis image
    colorbar()
    title('T1 map')

    figure(3001)
    imagesc(output_MRF_match_ssepg.T2_map)
    axis image
    colorbar()
    title('T2 map')

    figure(3002)
    imagesc(abs(output_MRF_match_ssepg.M0_map))
    axis image
    colorbar()
    title('Magnitude of M0')

    % save MRF match w/o fat sep

    save([dir_out fn_MRF_proc_no_fat_sep_epg],'output_MRF_match_epg')
    save([dir_out fn_MRF_proc_no_fat_sep_pepg],'output_MRF_match_pepg')
    save([dir_out fn_MRF_proc_no_fat_sep_ssepg],'output_MRF_match_ssepg')

end


