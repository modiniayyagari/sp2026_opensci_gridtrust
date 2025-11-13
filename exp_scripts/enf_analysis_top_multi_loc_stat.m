%% ENF Multi-Location Validation (Experiment 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   Quantifies cross-location separability of ENF signatures by comparing a sensed
%   ambient-EM ENF (near a DC-powered PCB/FPGA) against mains references collected
%   across multiple power-grid regions (e.g., US East/Texas/West 60 Hz, Germany 50 Hz).
%   Produces per-folder and summarized correlation heatmaps, descriptive statistics
%   on r and Fisher z', hypothesis tests (two-sample t on z'), and parametric EER
%   estimates under Gaussian modeling (with optional conservative variance adjustment).
%
% ANALYSIS OVERVIEW
% 1) INPUTS
%    • Root experiment directory containing trace-set folders {T01, T02, …}.
%    • Each trace-set folder holds ?2 time-aligned WAV files drawn from:
%        - fpga_em_trace_dc_*.wav      (ambient EM sensed near the FPGA board)
%        - mains_pow_trace_ac_*.wav    (grid mains references across locations)
%    • Filenames encode grid/region (egrid/tgrid/wgrid, DE_50, etc.) for grouping.
%
% 2) ENF EXTRACTION & CORRELATION
%    • STFT around nominal 50/60 Hz harmonics (configurable nfft/frame/overlap).
%    • Weighted-harmonic ENF estimator per file.
%    • Pearson correlation r between sensed ENF and reference ENF for every file pair.
%    • Aggregation of pairwise r into (a) per-folder matrices and (b) per-grid summaries.
%    • Fisher transform z = atanh(r) for parametric inference and modeling.
%
% 3) STATISTICS
%    • Descriptives on r and z (N, mean, SD, SEM, 95% CI) for intra-grid vs inter-grid.
%    • Effect sizes & separation diagnostics:
%        - Cohen’s d (on z), Cliff’s delta / overlap / margin (on r).
%    • Hypothesis test: two-sample t-test on z (intra-grid > inter-grid; Welch unequal-var).
%    • Parametric ROC point: Equal-Error Rate (EER) via Gaussian class model on z:
%        - ?* = (?0·?1 + ?1·?0) / (?0 + ?1), with FPR(?*) ? FNR(?*) = EER.
%        - Optional “conservative” model inflates ? for the genuine (intra) class
%          using temporal-reliability SDs from Experiment 3.
%    • Operating points at fixed errors (optional):
%        - Fix FPR = ? ? ? = ?0 + ?0·??¹(1??), then compute resulting FNR.
%        - Fix FNR = ? ? ? = ?1 + ?1·??¹(?), then compute resulting FPR.
%
% 4) VISUALIZATION & REPORTING
%    • Per-folder lower-triangle correlation heatmaps with annotated r.
%    • Per-file average heatmaps (within region) and grid-level summary matrices.
%    • Histograms/boxplots of r and z; optional z-PDF visualization with overlap/EER.
%    • Color map: blue = high correlation (intra-grid), orange = low correlation (inter-grid).
%
% 5) OUTPUTS
%    • Logs            ? exp_logs/
%    • Figures/Tables  ? exp_results/multiloc_validation_1/
%       - avg_file_heatmap_DE.svg / avg_file_heatmap_US.svg
%       - summary_correlation_heatmap_DE.svg / _US.svg
%       - r/z boxplots & histograms; z-distribution/EER plot (optional)
%    • Console/Log summaries of descriptives, tests, and EER/threshold metrics.
%
% Usage:
%   • Set config.baseDir to the MULTI root (contains T01/T02/…).
%   • Ensure proc_enf_analysis.m (or .p) and helpers are on the MATLAB path.
%   • Toggle outputs via config flags:
%       SCR_CFG_EN_CORR_VAL_MAT_PER_trace_set_FOLDER[_AVG], SCR_CFG_EN_CORR_VAL_MAT_PER_GRIDFREQ_AVG,
%       SCR_CFG_DISP_DESCRIPTIVE_STATS_RES, SCR_CFG_EN_TTEST.
%   • (Optional) Set conservative ? from Experiment 3 before EER computation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear; clc;

%% 1) Script Configuration
config = struct();

%Main Path Configuration
config.baseDir = '../exp_inputs/MULTI/';

%Feature Flags
config.SCR_CFG_DISP_RAW_PAIR_WISE_CORR_VALS             = false;
config.SCR_CFG_DISP_DESCRIPTIVE_STATS_RES               = true;
config.SCR_CFG_EN_CORR_VAL_MAT_PER_trace_set_FOLDER     = false;
config.SCR_CFG_EN_CORR_VAL_MAT_PER_trace_set_FOLDER_AVG = true;
config.SCR_CFG_EN_CORR_VAL_MAT_PER_GRIDFREQ_AVG         = true;
config.SCR_CFG_EN_TTEST                                 = true;

%Spectrogram / ENF extraction parameters
config.nominal_freq_arr = [50 60];
frame_size_arr          = (1:12)*1000;
config.frame_size       = frame_size_arr(8);
nfft_arr                = 2.^(10:20);
config.nfft             = nfft_arr(6);
overlap_size_arr        = 0:0.1:0.9;
config.overlap_size     = overlap_size_arr(1)*config.frame_size;

%% 2) Script Setup
[fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders();
config.results_dir = results_dir;

fprintf('Starting ENF Multi-Location Validation Analysis...\n');
fprintf(fileID, 'Starting ENF Multi-Location Validation Analysis... (Log file: %s)\n\n', log_filename);

%% 3) Main ENF Extraction and Analysis
fprintf('Finding all trace_set experiment folders...\n');
fprintf(fileID, 'Finding all trace_set experiment folders...\n');

trace_set_folders = find_trace_set_folders(config.baseDir, fileID);

fprintf('Found %d trace_set folders. Starting analysis...\n\n', length(trace_set_folders));
fprintf(fileID, 'Found %d trace_set folders. Starting analysis...\n\n', length(trace_set_folders));

[all_results, intra_grid_correlations, inter_grid_correlations] = ...
    run_pairwise_analysis(trace_set_folders, config, fileID);

%% 4) Experiment Trace Collection Summary
print_exp_analysis_summary(all_results, config, fileID);

%% 5) Statistical Analysis
% STEP 1A: Descriptive Results (r)
r_stats.intra = calculate_basic_stats(intra_grid_correlations);
r_stats.inter = calculate_basic_stats(inter_grid_correlations);

% STEP 1B: Fisher z'-transform and Descriptive Results
intra_grid_correlations_z = calculate_fisher_z_transform(intra_grid_correlations);
inter_grid_correlations_z = calculate_fisher_z_transform(inter_grid_correlations);
z_stats.intra = calculate_basic_stats(intra_grid_correlations_z);
z_stats.inter = calculate_basic_stats(inter_grid_correlations_z);

% STEP 1C: Comparative Statistics (Effect Sizes & Separation)
r_stats.cliff_delta = sub_calc_cliff_delta(intra_grid_correlations, inter_grid_correlations);
r_stats.overlap     = sub_calc_overlap(intra_grid_correlations, inter_grid_correlations);
r_stats.margin      = sub_calc_margin(intra_grid_correlations, inter_grid_correlations);
z_stats.cohens_d    = sub_calc_cohens_d(intra_grid_correlations_z, inter_grid_correlations_z);

if config.SCR_CFG_DISP_DESCRIPTIVE_STATS_RES
    fprintf('\n\n--- STEP 1A: Descriptive Results ---\n');
    fprintf(fileID, '\n\n--- STEP 1A: Descriptive Results ---\n');   
    
    fprintf('\n\n--- STEP 1A: Basic Statistical Values for Raw Intra and Inter-Grid Correlation Values ---\n');
    fprintf(fileID, '\n\n--- STEP 1A: Basic Statistical Values for Raw Intra and Inter-Grid Correlation Values ---\n');    
    print_stats_table(r_stats.intra, r_stats.inter, fileID, 'Intra-Grid (r)', 'Inter-Grid (r)');

    plot_descriptive_stats_results_plot_boxplot( ...
        intra_grid_correlations, inter_grid_correlations, r_stats.intra, config.results_dir, fileID, 'r');   
    plot_descriptive_stats_results_plot_histogram( ...
        intra_grid_correlations, inter_grid_correlations, r_stats.intra, config.results_dir, fileID, 'r');

    fprintf('\n\n--- STEP 1B: Basic Statistical Values for Z-transformed Intra and Inter-Grid Correlation Values ---\n');
    fprintf(fileID, '\n\n--- STEP 1B: Basic Statistical Values for Z-transformed Intra and Inter-Grid Correlation Values ---\n');
    print_stats_table(z_stats.intra, z_stats.inter, fileID, 'Intra-Grid (z'')', 'Inter-Grid (z'')');

    plot_descriptive_stats_results_plot_boxplot( ...
        intra_grid_correlations_z, inter_grid_correlations_z, z_stats.intra, config.results_dir, fileID, 'z');
    plot_descriptive_stats_results_plot_histogram( ...
        intra_grid_correlations_z, inter_grid_correlations_z, z_stats.intra, config.results_dir, fileID, 'z');

    fprintf('\nGenerating and saving z-score distribution plot...\n');
    fprintf(fileID, '\nGenerating and saving z-score distribution plot...\n');
    try
        plot_z_score_distributions(z_stats, config.results_dir, fileID);
    catch ME_z_plot
        fprintf('    Error generating z-score distribution plot: %s\n', ME_z_plot.message);
        fprintf(fileID, '    Error generating z-score distribution plot: %s\n', ME_z_plot.message);
    end

    fprintf('\n\n--- STEP 1C: Comparison Metrics (Effect Size & Separation) ---\n');
    fprintf(fileID, '\n\n--- STEP 1C: Comparison Metrics (Effect Size & Separation) ---\n');    
    fprintf(        '%-25s | %18s\n', 'Metric', 'Value');
    fprintf(fileID, '%-25s | %18s\n', 'Metric', 'Value');
    fprintf(        '%-25s | %18s\n', '-------------------------', '------------------');
    fprintf(fileID, '%-25s | %18s\n', '-------------------------', '------------------');    
    fprintf(        '%-25s | %18.4f\n', 'Cohen''s d (from z'')', z_stats.cohens_d);
    fprintf(fileID, '%-25s | %18.4f\n', 'Cohen''s d (from z'')', z_stats.cohens_d);    
    fprintf(        '%-25s | %18.4f\n', 'Cliff''s delta (from r)', r_stats.cliff_delta);
    fprintf(fileID, '%-25s | %18.4f\n', 'Cliff''s delta (from r)', r_stats.cliff_delta);    
    fprintf(        '%-25s | %18.4f\n', 'Overlap Coefficient (r)', r_stats.overlap);
    fprintf(fileID, '%-25s | %18.4f\n', 'Overlap Coefficient (r)', r_stats.overlap);    
    fprintf(        '%-25s | %18.4f\n', 'Separation Margin (r)', r_stats.margin);
    fprintf(fileID, '%-25s | %18.4f\n', 'Separation Margin (r)', r_stats.margin);
        
    fprintf('\n\n--- STEP 1 Complete ---\n');
    fprintf(fileID, '\n\n--- STEP 1 Complete ---\n');
else
    fprintf('\n\n--- STEP 1: Descriptive Results (Skipped) ---\n');
    fprintf(fileID, '\n\n--- STEP 1: Descriptive Results (Skipped) ---\n');
end

% STEP 2A: Per-File Average Summary Heatmap
if config.SCR_CFG_EN_CORR_VAL_MAT_PER_trace_set_FOLDER_AVG
    generate_averaged_trace_set_heatmap(all_results, config.results_dir, fileID);
else
    fprintf('\n\nSTEP 2A: Per-File Average Summary Heatmap (Skipped)\n');
    fprintf(fileID, '\n\nSTEP 2A: Per-File Average Summary Heatmap (Skipped)\n');
end

% STEP 2B: Summarized Grid Heatmap
if config.SCR_CFG_EN_CORR_VAL_MAT_PER_GRIDFREQ_AVG
    generate_averaged_grid_freq_heatmap(all_results, config.results_dir, fileID);
else
    fprintf('\n\nSTEP 2B: Summarized Correlation Matrix Heatmap (Skipped)\n');
    fprintf(fileID, '\n\nSTEP 2B: Summarized Correlation Matrix Heatmap (Skipped)\n');
end

% STEP 3: Formal Hypothesis Testing (t-test)
if config.SCR_CFG_EN_TTEST
    enable_t_test_analysis(intra_grid_correlations_z, inter_grid_correlations_z, fileID);
else
    fprintf('\n\nSTEP 3: Formal Hypothesis Testing (Skipped)\n');
    fprintf(fileID, '\n\nSTEP 3: Formal Hypothesis Testing (Skipped)\n');
end

% STEP 4: Parametric EER (Gaussian approximation on Fisher z')
eer_stats = compute_parametric_eer_gaussian(z_stats, fileID);

fprintf('\n\nAll analysis complete. Log file saved.\n');
fprintf(fileID, '\n\nAll analysis complete.\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedures for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders()
    log_dir = 'exp_logs';
    results_dir = 'exp_results/enf_analysis_top_multi_loc_stat';
    
    if ~exist(log_dir, 'dir'), mkdir(log_dir); end
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    
    [~, script_name] = fileparts(mfilename('fullpath'));
    if isempty(script_name)
        script_name = 'enf_analysis_top_multi_loc_stat';
    end
    log_filename = fullfile(log_dir, [script_name, '_log.txt']);
    
    if exist(log_filename, 'file'), delete(log_filename); end
    
    fileID = fopen(log_filename, 'w');
    if fileID == -1
        error('Could not open log file %s for writing.', log_filename);
    end
    
    cleanupObj = onCleanup(@() fclose(fileID));
end

function cmap = setup_enf_colormap()
% ENF diverging colormap tuned for correlation heatmaps.
% [-0.060, 0.40] -> orange tones
%  0.50         -> white (pivot)
% [0.60, 1.00]  -> blue tones
% Extra contrast & resolution in [0.90, 1.00] with gamma shaping.

    vmin        = -0.10;
    vwhite      =  0.50;
    vblue_start =  0.60;   %#ok<NASGU> (kept for clarity)
    vfocus      =  0.99;   % start of high-resolution blue zone
    vmax        =  1.00;

    % Colors
    orange     = [253,174, 97]/255;   % warm orange for low values
    white      = [1,1,1];            % neutral center
    blue_mid   = [144,202,249]/255;  % light/medium blue
    blue_dark  = [66,133,244]/255;  % lighter "deep" blue so black text is legible

    % High-res samples for [0.90, 1.00]
    hi_step = 0.001;
    n_hi    = ceil((vmax - vfocus)/hi_step) + 1;

    % Total samples in colormap
    N_total = 2048 + n_hi;

    % Allocate between segments [vmin,vwhite], [vwhite,vfocus]
    base_span = (vwhite - vmin) + (vfocus - vwhite); % = 1.0 here
    N_base    = N_total - n_hi;

    n_low = max(80, round(N_base * (vwhite - vmin)/base_span)); % orange -> white
    n_mid = max(80, N_base - n_low);                            % white -> blue_mid

    % Build value grid
    vals_low  = linspace(vmin,   vwhite, n_low);                  % includes white
    vals_mid  = linspace(vwhite, vfocus, n_mid+1); vals_mid(1) = [];
    vals_high = linspace(vfocus, vmax,   n_hi+1);  vals_high(1) = [];

    vals = [vals_low, vals_mid, vals_high];

    % Gamma shaping for [vfocus, vmax] to emphasize differences near 1.0
    gamma_blue = 3;

    cmap = zeros(numel(vals), 3);

    for k = 1:numel(vals)
        v = vals(k);

        if v <= vwhite
            % Orange -> White
            t = (v - vmin) / (vwhite - vmin);
            t = max(min(t,1),0);
            cmap(k,:) = orange + t*(white - orange);

        elseif v <= vfocus
            % White -> Medium Blue
            t = (v - vwhite) / (vfocus - vwhite);
            t = max(min(t,1),0);
            cmap(k,:) = white + t*(blue_mid - white);

        else
            % Medium Blue -> Lighter Dark Blue with gamma shaping
            t_lin = (v - vfocus) / (vmax - vfocus);
            t = max(min(t_lin,1),0).^gamma_blue;
            cmap(k,:) = blue_mid + t*(blue_dark - blue_mid);
        end
    end
end

function trace_set_folders = find_trace_set_folders(baseDir, fileID)
    all_paths_str = genpath(baseDir);
    all_paths = strsplit(all_paths_str, pathsep);
    
    trace_set_folders_list = struct('name', {}, 'folder', {}, 'isdir', {});
    
    for i = 1:length(all_paths)
        current_path = all_paths{i};
        if isempty(current_path), continue; end
        
        [parent_folder, folder_name] = fileparts(current_path);
        
        if ~isempty(folder_name) && startsWith(folder_name, 'T')
            s.name = folder_name;
            s.folder = parent_folder;
            s.isdir = true;
            trace_set_folders_list(end+1) = s;
        end
    end
    trace_set_folders = trace_set_folders_list;
end

function [all_results, intra_grid_correlations, inter_grid_correlations] = ...
    run_pairwise_analysis(trace_set_folders, config, fileID)

    intra_grid_correlations = [];
    inter_grid_correlations = [];
    all_results = struct('Folder', {}, 'File_A', {}, 'Grid_A', {}, ...
                         'File_B', {}, 'Grid_B', {}, 'Correlation', {});

    for k = 1:length(trace_set_folders)
        trace_set_path = fullfile(trace_set_folders(k).folder, trace_set_folders(k).name);    
        trace_set_path = strrep(trace_set_path, '\', '/');
        relative_baseDir = strrep(trace_set_path, config.baseDir, '');    
        safe_folder_prefix = strrep(relative_baseDir, '/', '_');    
        if startsWith(safe_folder_prefix, '_'), safe_folder_prefix = safe_folder_prefix(2:end); end

        fprintf('Processing Folder: %s\n', trace_set_path);
        fprintf(fileID, 'Processing Folder: %s\n', trace_set_path);
        
        wav_files = dir(fullfile(trace_set_path, '*.wav'));
        num_files = length(wav_files);
        
        if num_files < 2
            fprintf('  Skipping (found %d .wav files, need at least 2 for pairs).\n\n', num_files);
            fprintf(fileID, '  Skipping (found %d .wav files, need at least 2 for pairs).\n\n', num_files);
            continue;
        end
        
        fprintf('  Found %d traces. Running %d pairwise correlations...\n', num_files, nchoosek(num_files, 2));
        fprintf(fileID, '  Found %d traces. Running %d pairwise correlations...\n', num_files, nchoosek(num_files, 2));

        for i = 1:num_files
            for j = (i + 1):num_files
                
                file_A_name = wav_files(i).name;
                file_B_name = wav_files(j).name;
                
                full_path_A = fullfile(trace_set_path, file_A_name);
                full_path_B = fullfile(trace_set_path, file_B_name);
                
                id_A = get_file_identity(file_A_name, trace_set_path, config.nominal_freq_arr);
                id_B = get_file_identity(file_B_name, trace_set_path, config.nominal_freq_arr);

                corr_val = proc_enf_analysis(full_path_A, full_path_B, ...
                                             config.nfft, config.frame_size, config.overlap_size, ...
                                             id_A.params.harmonics_arr, id_A.params.nominal_freq, ...
                                             id_B.params.harmonics_arr, id_B.params.nominal_freq, ...
                                             id_A.params.est_method, id_A.params.est_freq, id_A.params.est_spec_comb_harmonics, ...
                                             id_B.params.est_method, id_B.params.est_freq, id_B.params.est_spec_comb_harmonics, ...
                                             id_A.params.plot_title, id_B.params.plot_title, false);
                                         
                
                new_result.Folder      = safe_folder_prefix;
                new_result.File_A      = file_A_name;
                new_result.Grid_A      = id_A.grid;
                new_result.File_B      = file_B_name;
                new_result.Grid_B      = id_B.grid;
                new_result.Correlation = corr_val;
                all_results(end+1) = new_result; %#ok<AGROW>
                
                if strcmp(id_A.grid, id_B.grid)
                    intra_grid_correlations(end+1) = corr_val; %#ok<AGROW>
                else
                    inter_grid_correlations(end+1) = corr_val; %#ok<AGROW>
                end
                
                log_msg = sprintf('    > Compared (%s) vs. (%s): r = %.4f [Group: %s]\n', ...
                                  file_A_name, file_B_name, corr_val, ...
                                  iif(strcmp(id_A.grid, id_B.grid), 'Intra-Grid', 'Inter-Grid'));
                fprintf('%s', log_msg);
                fprintf(fileID, '%s', log_msg);
            end
        end
        
        if config.SCR_CFG_EN_CORR_VAL_MAT_PER_trace_set_FOLDER
            generate_per_trace_set_heatmap(wav_files, all_results, safe_folder_prefix, ...
                                           config.results_dir, fileID);        
        end
        
        fprintf('  Finished folder.\n\n');
        fprintf(fileID, '  Finished folder.\n\n');
    end

    fprintf('Analysis Complete\n\n');
    fprintf(fileID, 'Analysis Complete\n\n');
end

function print_exp_analysis_summary(all_results, config, fileID)

    fprintf('\n\n Full Experiment Summary \n');
    fprintf(fileID, '\n\n Full Experiment Summary \n');

    try
        de_results_idx = ~cellfun('isempty', strfind({all_results.Folder}, 'DE_'));
        de_results = all_results(de_results_idx);
        us_results_idx = ~cellfun('isempty', strfind({all_results.Folder}, 'US_'));
        us_results = all_results(us_results_idx);

        all_folders = {all_results.Folder};
        de_folders  = unique(all_folders(de_results_idx));
        us_folders  = unique(all_folders(us_results_idx));

        fprintf('  Total T0* folders analyzed: %d (DE_50: %d, US_60: %d)\n', ...
            length(de_folders) + length(us_folders), length(de_folders), length(us_folders));
        fprintf(fileID, '  Total T0* folders analyzed: %d (DE_50: %d, US_60: %d)\n', ...
            length(de_folders) + length(us_folders), length(de_folders), length(us_folders));

        count_pairs = @(results, fileA, fileB) sum( ...
            (strcmp({results.File_A}, fileA) & strcmp({results.File_B}, fileB)) | ...
            (strcmp({results.File_B}, fileA) & strcmp({results.File_A}, fileB)) );

        fprintf('\n  DE_50 Summary\n');
        fprintf(fileID, '\n  DE_50 Summary\n');

        f_de_em      = 'fpga_em_trace_dc_citya_lab.wav';
        f_de_dresden = 'mains_pow_trace_ac_dresden.wav';
        f_de_control = 'mains_pow_trace_ac_citya_lab.wav';
        
        count_de_1 = count_pairs(de_results, f_de_em, f_de_dresden);
        count_de_2 = count_pairs(de_results, f_de_em, f_de_control);

        fprintf('    Total pairs (fpga_citya_lab vs mains_dresden): %d\n', count_de_1);
        fprintf(fileID, '    Total pairs (fpga_citya_lab vs mains_dresden): %d\n', count_de_1);
        fprintf('    Total pairs (fpga_citya_lab vs mains_citya_lab): %d\n', count_de_2);
        fprintf(fileID, '    Total pairs (fpga_citya_lab vs mains_citya_lab): %d\n', count_de_2);

        de_intra_count = sum(strcmp({de_results.Grid_A}, {de_results.Grid_B}));
        de_inter_count = sum(~strcmp({de_results.Grid_A}, {de_results.Grid_B}));
        fprintf('    Total Intra-Grid Comparisons: %d\n', de_intra_count);
        fprintf(fileID, '    Total Intra-Grid Comparisons: %d\n', de_intra_count);
        fprintf('    Total Inter-Grid Comparisons: %d\n', de_inter_count);
        fprintf(fileID, '    Total Inter-Grid Comparisons: %d\n', de_inter_count);

        fprintf('\n  US_60 Summary\n');
        fprintf(fileID, '\n  US_60 Summary\n');

        f_us_em_lab = 'fpga_em_trace_dc_egrid_citya_lab.wav';
        f_us_ac_lab = 'mains_pow_trace_ac_egrid_citya_lab.wav';
        f_us_ac_home = 'mains_pow_trace_ac_egrid_citya_home.wav';
        f_us_ac_worc = 'mains_pow_trace_ac_egrid_worcester.wav';
        f_us_ac_rich = 'mains_pow_trace_ac_tgrid_richardson.wav';
        f_us_ac_tucs = 'mains_pow_trace_ac_wgrid_tucson.wav';

        us_pairs_to_count = {
            f_us_ac_lab,  'fpga_citya_lab vs mains_citya_lab';
            f_us_ac_home, 'fpga_citya_lab vs mains_citya_home';
            f_us_ac_worc, 'fpga_citya_lab vs mains_worcester';
            f_us_ac_rich, 'fpga_citya_lab vs mains_richardson';
            f_us_ac_tucs, 'fpga_citya_lab vs mains_tucson'
        };
    
        for i = 1:size(us_pairs_to_count, 1)
            fileB = us_pairs_to_count{i, 1};
            label = us_pairs_to_count{i, 2};
            count = count_pairs(us_results, f_us_em_lab, fileB);
            fprintf('    Total pairs (%s): %d\n', label, count);
            fprintf(fileID, '    Total pairs (%s): %d\n', label, count);
        end

        us_intra_count = sum(strcmp({us_results.Grid_A}, {us_results.Grid_B}));
        us_inter_count = sum(~strcmp({us_results.Grid_A}, {us_results.Grid_B}));
        fprintf('    Total Intra-Grid Comparisons: %d\n', us_intra_count);
        fprintf(fileID, '    Total Intra-Grid Comparisons: %d\n', us_intra_count);
        fprintf('    Total Inter-Grid Comparisons: %d\n', us_inter_count);
        fprintf(fileID, '    Total Inter-Grid Comparisons: %d\n', us_inter_count);
    
    catch ME_summary
        fprintf('  Error generating full experiment summary: %s\n', ME_summary.message);
        fprintf(fileID, '  Error generating full experiment summary: %s\n', ME_summary.message);
    end
    
    if config.SCR_CFG_DISP_RAW_PAIR_WISE_CORR_VALS
        fprintf('\n\nFull Results Log\n');
        fprintf(fileID, '\n\nFull Results Log\n');
        T_str = evalc('disp(struct2table(all_results))');
        fprintf('%s\n', T_str);
        fprintf(fileID, '%s\n', T_str);
    end
    
    fprintf('\n\n ################################################## \n');
    fprintf(fileID, '\n\n ################################################## \n');
end

function z_transform_data = calculate_fisher_z_transform(r_data)
    r_bounded = r_data;
    r_bounded(r_bounded >= 1.0)  = 0.9999999;
    r_bounded(r_bounded <= -1.0) = -0.9999999;
    z_transform_data = atanh(r_bounded);
end

function stats = calculate_basic_stats(data_array)
    if ~isempty(data_array)
        stats.N    = length(data_array);
        stats.mean = mean(data_array);
        stats.median = median(data_array);
        stats.std  = std(data_array);
        stats.min  = min(data_array);
        stats.max  = max(data_array);
        stats.var  = var(data_array);
        stats.skew = skewness(data_array);
        stats.kurt = kurtosis(data_array);
        stats.sem  = stats.std / sqrt(stats.N);
        stats.ci95_margin = 1.96 * stats.sem;
    else
        stats.N = 0;
        stats.mean = NaN; stats.median = NaN; stats.std = NaN;
        stats.min = NaN; stats.max = NaN; stats.var = NaN;
        stats.skew = NaN; stats.kurt = NaN; stats.sem = NaN;
        stats.ci95_margin = NaN;
    end
end

function cliff_d = sub_calc_cliff_delta(g1, g2)
    if isempty(g1) || isempty(g2)
        cliff_d = NaN;
        return;
    end
    n1 = length(g1);
    n2 = length(g2);
    g1_matrix = repmat(g1(:), 1, n2);
    g2_matrix = repmat(g2(:)', n1, 1);
    comparisons = sign(g1_matrix - g2_matrix);
    cliff_d = sum(comparisons(:)) / (n1 * n2);
end

function ovl = sub_calc_overlap(g1, g2)
    if isempty(g1) || isempty(g2)
        ovl = NaN;
        return;
    end
    try
        min_val = min(min(g1), min(g2));
        max_val = max(max(g1), max(g2));
        x_range = linspace(min_val, max_val, 1000);
        
        [f1, ~] = ksdensity(g1, x_range);
        [f2, ~] = ksdensity(g2, x_range);
        
        f1 = f1 / trapz(x_range, f1);
        f2 = f2 / trapz(x_range, f2);
        
        overlap_curve = min(f1, f2);
        ovl = trapz(x_range, overlap_curve);
    catch
        ovl = NaN;
    end
end

function margin = sub_calc_margin(g1, g2)
    if isempty(g1) || isempty(g2)
        margin = NaN;
        return;
    end
    margin = min(g1) - max(g2);
end

function cohens_d = sub_calc_cohens_d(g1, g2)
    if isempty(g1) || isempty(g2)
        cohens_d = NaN;
        return;
    end
    n1 = length(g1);
    n2 = length(g2);
    mean1 = mean(g1);
    mean2 = mean(g2);
    std1 = std(g1);
    std2 = std(g2);
    
    s_pooled = sqrt(((n1-1)*std1^2 + (n2-1)*std2^2) / (n1 + n2 - 2));
    
    if s_pooled == 0
        cohens_d = Inf;
    else
        cohens_d = (mean1 - mean2) / s_pooled;
    end
end

function print_stats_table(intra_stats, inter_stats, fileID, intra_label, inter_label)
    fprintf(        '%-12s | %18s | %18s\n', 'Statistic', intra_label, inter_label);
    fprintf(fileID, '%-12s | %18s | %18s\n', 'Statistic', intra_label, inter_label);
    
    fprintf(        '%-12s | %18d | %18d\n', 'Count (N)', intra_stats.N, inter_stats.N);
    fprintf(fileID, '%-12s | %18d | %18d\n', 'Count (N)', intra_stats.N, inter_stats.N);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Mean', intra_stats.mean, inter_stats.mean);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Mean', intra_stats.mean, inter_stats.mean);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Median', intra_stats.median, inter_stats.median);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Median', intra_stats.median, inter_stats.median);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Std. Dev.', intra_stats.std, inter_stats.std);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Std. Dev.', intra_stats.std, inter_stats.std);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Variance', intra_stats.var, inter_stats.var);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Variance', intra_stats.var, inter_stats.var);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Skewness', intra_stats.skew, inter_stats.skew);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Skewness', intra_stats.skew, inter_stats.skew);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Kurtosis', intra_stats.kurt, inter_stats.kurt);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Kurtosis', intra_stats.kurt, inter_stats.kurt);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Std. Error', intra_stats.sem, inter_stats.sem);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Std. Error', intra_stats.sem, inter_stats.sem);
    
    intra_ci_str = sprintf('[%.4f, %.4f]', intra_stats.mean - intra_stats.ci95_margin, intra_stats.mean + intra_stats.ci95_margin);
    inter_ci_str = sprintf('[%.4f, %.4f]', inter_stats.mean - inter_stats.ci95_margin, inter_stats.mean + inter_stats.ci95_margin);
    
    fprintf(        '%-12s | %18s | %18s\n', '95% CI', intra_ci_str, inter_ci_str);
    fprintf(fileID, '%-12s | %18s | %18s\n', '95% CI', intra_ci_str, inter_ci_str);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Min', intra_stats.min, inter_stats.min);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Min', intra_stats.min, inter_stats.min);
    
    fprintf(        '%-12s | %18.4f | %18.4f\n', 'Max', intra_stats.max, inter_stats.max);
    fprintf(fileID, '%-12s | %18.4f | %18.4f\n', 'Max', intra_stats.max, inter_stats.max);
end

function plot_descriptive_stats_results_plot_boxplot(intra_grid_data, inter_grid_data, intra_stats, results_dir, fileID, value_type)

    if strcmp(value_type, 'z')
        plot_name   = 'Z-Transformed Correlation Box Plot';
        y_label     = 'Fisher''s z''-transform';
        title_suffix = '(z''-values)';
        file_suffix = '_z_transform';
    else
        plot_name   = 'Correlation Box Plot';
        y_label     = 'Pearson Correlation (r)';
        title_suffix = '(r-values)';
        file_suffix = '_r_value';
    end

    try
        fig_box = figure('Name', plot_name, 'Visible', 'off', 'Position', [100, 100, 1000, 600]);

        h_main = subplot(1, 2, 1);
        
        data = [intra_grid_data'; inter_grid_data']; 
        groups = [repmat({'Intra-Grid'}, length(intra_grid_data), 1); ...
                  repmat({'Inter-Grid'}, length(inter_grid_data), 1)]; 
        
        boxplot(h_main, data, groups, 'Notch', 'on');
        title(h_main, ['Comparison of Intra-Grid vs. Inter-Grid Correlations ' title_suffix]); 
        ylabel(h_main, y_label); 
        grid(h_main, 'on');
        
        h_zoom = subplot(1, 2, 2);
        
        if ~isempty(intra_grid_data) 
            boxplot(h_zoom, intra_grid_data', 'Notch', 'on', 'Labels', {'Intra-Grid (Zoomed)'}, ...
                    'Colors', [144, 202, 249]/255); 
            
            margin = (intra_stats.max - intra_stats.min) * 0.1;
            if margin == 0 || isnan(margin)
                margin = 0.001;
            end
            
            ylim(h_zoom, [intra_stats.min - margin, intra_stats.max + margin]);
            
            title(h_zoom, ['Zoomed: Intra-Grid Correlations ' title_suffix]); 
            ylabel(h_zoom, y_label); 
            grid(h_zoom, 'on');
        else
            title(h_zoom, 'Zoomed: Intra-Grid Correlations (No Data)');
            set(h_zoom, 'XTick', [], 'YTick', []);
        end
        
        save_filename = ['intra_vs_inter_grid_corr_vals_boxplot' file_suffix];
        saveas(fig_box, fullfile(results_dir, [save_filename '.svg']));
        fprintf('\nSaved %s.svg\n', save_filename);
        fprintf(fileID, '\nSaved %s.svg\n', save_filename);
        close(fig_box);
    catch ME
        fprintf('Error generating box plot (%s): %s\n', value_type, ME.message);
        fprintf(fileID, 'Error generating box plot (%s): %s\n', value_type, ME.message);
    end
end

function plot_descriptive_stats_results_plot_histogram(intra_grid_data, inter_grid_data, intra_stats, results_dir, fileID, value_type)

    if strcmp(value_type, 'z')
        plot_name    = 'Z-Transformed Correlation Histograms';
        x_label      = 'Fisher''s z''-transform';
        title_suffix = '(z''-values)';
        file_suffix  = '_z_transform';
        inter_xlim   = 'auto';
    else
        plot_name    = 'Correlation Histograms';
        x_label      = 'Pearson Correlation (r)';
        title_suffix = '(r-values)';
        file_suffix  = '_r_value';
        inter_xlim   = [-1, 1];
    end

    try
        fig_hist = figure('Name', plot_name, 'Visible', 'off', 'Position', [100, 100, 1000, 500]);
        
        h_inter = subplot(1, 2, 1);
        if ~isempty(inter_grid_data) 
            histogram(h_inter, inter_grid_data, 'BinMethod', 'auto', ...
                      'FaceColor', [229, 115, 115]/255, 'EdgeColor', 'w'); 
        end
        title(h_inter, ['Inter-Grid Correlations ' title_suffix]); 
        xlabel(h_inter, x_label); 
        ylabel(h_inter, 'Count');
        if ~ischar(inter_xlim), xlim(h_inter, inter_xlim); end
        grid(h_inter, 'on');
        
        h_intra = subplot(1, 2, 2);
        if ~isempty(intra_grid_data) 
            histogram(h_intra, intra_grid_data, 'BinMethod', 'auto', ...
                      'FaceColor', [144, 202, 249]/255, 'EdgeColor', 'w'); 
            
            margin = (intra_stats.max - intra_stats.min) * 0.1;
            if margin == 0 || isnan(margin)
                margin = 0.001;
            end
            
            xlim(h_intra, [intra_stats.min - margin, intra_stats.max + margin]);
        end
        title(h_intra, ['Intra-Grid Correlations (Zoomed, ' title_suffix ')']); 
        xlabel(h_intra, x_label); 
        ylabel(h_intra, 'Count');
        grid(h_intra, 'on');        
        
        save_filename = ['intra_vs_inter_grid_corr_vals_histogram' file_suffix];
        saveas(fig_hist, fullfile(results_dir, [save_filename '.svg']));
        fprintf('\nSaved %s.svg\n', save_filename);
        fprintf(fileID, '\nSaved %s.svg\n', save_filename);
        close(fig_hist);
    catch ME
        fprintf('Error generating histograms (%s): %s\n', value_type, ME.message);
        fprintf(fileID, 'Error generating histograms (%s): %s\n', value_type, ME.message);
    end
end

function generate_per_trace_set_heatmap(wav_files, all_results, safe_folder_prefix, results_dir, fileID)

    try
        fprintf('  Generating per-folder heatmap...\n');
        fprintf(fileID, '  Generating per-folder heatmap...\n');
        
        folder_files = {wav_files.name};
        n = numel(folder_files);
        labels = cell(1,n);
        for i = 1:n
            labels{i} = get_corrmat_label(folder_files{i});
        end
        
        % Build lower-triangular correlation matrix
        M = nan(n);
        this_idx = strcmp({all_results.Folder}, safe_folder_prefix);
        folder_results = all_results(this_idx);
        for i = 1:n
            M(i,i) = 1.0;
            for j = 1:(i-1)
                fa = folder_files{i};
                fb = folder_files{j};
                idx = find(strcmp({folder_results.File_A},fa) & strcmp({folder_results.File_B},fb),1);
                if isempty(idx)
                    idx = find(strcmp({folder_results.File_A},fb) & strcmp({folder_results.File_B},fa),1);
                end
                if ~isempty(idx)
                    M(i,j) = folder_results(idx).Correlation;
                end
            end
        end
        
        % True value range for mapping
        vmin = -0.06;
        vmax =  1.0;
        gamma = 0.35;        % < 1 => more contrast near 1
        r_center = (vmin + vmax)/2;
        
        % Map real r -> [0,1] -> gamma -> colormap index
        T  = (M - vmin) / (vmax - vmin);
        T  = max(min(T,1),0);
        Tg = T .^ gamma;
        
        % Center of diverging colormap in [0,1] after gamma
        t_center = (r_center - vmin) / (vmax - vmin);
        t_center_g = t_center ^ gamma;
        
        numColors = 2048;
        cmap = setup_diverging_colormap(0, 1, t_center_g, numColors);
        
        fig = figure('Name',['Heatmap ' safe_folder_prefix], ...
                     'Visible','off', 'Position',[100 100 800 650]);
        
        imagesc(Tg, 'AlphaData', ~isnan(M));
        set(gca,'Color',[1 1 1]);
        axis equal tight;
        colormap(cmap);
        caxis([0 1]);
        
        % Colorbar with ticks relabeled to true r
        cb = colorbar;
        nt = 6;
        ticks = linspace(0,1,nt);
        cb.Ticks = ticks;
        true_vals = vmin + (ticks .^ (1/gamma)) * (vmax - vmin);
        cb.TickLabels = arrayfun(@(x) sprintf('%.2f',x), true_vals, 'UniformOutput',false);
        ylabel(cb,'Correlation');
        
        xticks(1:n); yticks(1:n);
        xticklabels(labels); yticklabels(labels);
        xtickangle(45);
        
        for r = 1:n
            for c = 1:n
                if ~isnan(M(r,c))
                    text(c, r, sprintf('%.3f', M(r,c)), ...
                        'HorizontalAlignment','center', ...
                        'FontSize',9, 'Color','k');
                end
            end
        end
        
        box on;
        title(['Per-Folder Correlation: ' strrep(safe_folder_prefix,'_',' ')]);
        xlabel('Trace Type');
        ylabel('Trace Type');
        
        caption_str = generate_heatmap_caption(labels);
        if ~isempty(caption_str)
            add_heatmap_caption(fig, caption_str);
        end
        
        fname = fullfile(results_dir, ['corrmat_' safe_folder_prefix '.svg']);
        saveas(fig, fname);
        fprintf('  Saved %s\n', fname);
        fprintf(fileID, '  Saved %s\n', fname);
        close(fig);
        
    catch ME
        fprintf('  Error generating per-folder heatmap: %s\n', ME.message);
        fprintf(fileID, '  Error generating per-folder heatmap: %s\n', ME.message);
    end
end

function generate_averaged_trace_set_heatmap(all_results, results_dir, fileID)

    fprintf('\n\nSTEP 2A: Per-File Average Summary Heatmap\n');
    fprintf(fileID, '\n\nSTEP 2A: Per-File Average Summary Heatmap\n');

    vmin = -0.060;
    vmax =  1.00;
    cmap = setup_enf_colormap();

    %% DE 50Hz
    try
        fprintf('  Generating DE_50 Average Per-File Heatmap...\n');
        fprintf(fileID, '  Generating DE_50 Average Per-File Heatmap...\n');

        de_idx = ~cellfun('isempty', strfind({all_results.Folder}, 'DE_'));
        de_results = all_results(de_idx);

        if ~isempty(de_results)
            de_raw = unique([{de_results.File_A}, {de_results.File_B}]);
            n = numel(de_raw);
            labels = cell(1,n);
            for i=1:n, labels{i} = get_corrmat_label(de_raw{i}); end

            M = nan(n);
            for i = 1:n
                M(i,i) = 1;
                for j = 1:(i-1)
                    a = de_raw{i}; b = de_raw{j};
                    idx1 = find(strcmp({de_results.File_A},a)&strcmp({de_results.File_B},b));
                    idx2 = find(strcmp({de_results.File_A},b)&strcmp({de_results.File_B},a));
                    idx  = [idx1 idx2];
                    if ~isempty(idx), M(i,j) = mean([de_results(idx).Correlation]); end
                end
            end

            fig = figure('Name','DE Avg Per-File Heatmap','Visible','off');
            imagesc(M, 'AlphaData', ~isnan(M));
            set(gca,'Color',[1 1 1]);
            axis equal tight;
            colormap(cmap);
            caxis([vmin vmax]);
            cb = colorbar; ylabel(cb,'Correlation');

            xticks(1:n); yticks(1:n);
            xticklabels(labels); yticklabels(labels);
            xtickangle(45);
            for r=1:n, for c=1:n
                if ~isnan(M(r,c))
                    text(c,r,sprintf('%.3f',M(r,c)), ...
                        'HorizontalAlignment','center','FontSize',9,'Color','k');
                end
            end, end
            box on;
            title('Average Per-File Correlation (DE 50Hz Experiment)');
            xlabel('Trace Type'); ylabel('Trace Type');

            saveas(fig, fullfile(results_dir,'avg_file_heatmap_DE.svg'));
            fprintf('  Saved avg_file_heatmap_DE.svg\n');
            fprintf(fileID,'  Saved avg_file_heatmap_DE.svg\n');
            close(fig);
        else
            fprintf('  No DE data found, skipping.\n');
            fprintf(fileID,'  No DE data found, skipping.\n');
        end
    catch ME
        fprintf('  Error generating DE average per-file heatmap: %s\n', ME.message);
        fprintf(fileID,'  Error generating DE average per-file heatmap: %s\n', ME.message);
    end

    %% US 60Hz (same style)
    try
        fprintf('  Generating US_60 Average Per-File Heatmap...\n');
        fprintf(fileID, '  Generating US_60 Average Per-File Heatmap...\n');

        us_idx = ~cellfun('isempty', strfind({all_results.Folder}, 'US_'));
        us_results = all_results(us_idx);

        if ~isempty(us_results)
            % desired label order
            pretty = {'EM','AR1','AR2','AR3','AR4','AR5'};
            all_raw = unique([{us_results.File_A}, {us_results.File_B}]);

            raw_ordered = cell(size(pretty));
            for i=1:numel(pretty)
                for j=1:numel(all_raw)
                    if strcmp(get_corrmat_label(all_raw{j}), pretty{i})
                        raw_ordered{i} = all_raw{j};
                        break;
                    end
                end
            end
            valid = ~cellfun('isempty', raw_ordered);
            raw = raw_ordered(valid);
            labels = pretty(valid);

            n = numel(labels);
            M = nan(n);
            for i=1:n
                M(i,i) = 1;
                for j=1:(i-1)
                    a = raw{i}; b = raw{j};
                    idx1 = find(strcmp({us_results.File_A},a)&strcmp({us_results.File_B},b));
                    idx2 = find(strcmp({us_results.File_A},b)&strcmp({us_results.File_B},a));
                    idx = [idx1 idx2];
                    if ~isempty(idx), M(i,j) = mean([us_results(idx).Correlation]); end
                end
            end

            fig = figure('Name','US Avg Per-File Heatmap','Visible','off');
            imagesc(M, 'AlphaData', ~isnan(M));
            set(gca,'Color',[1 1 1]);
            axis equal tight;
            colormap(cmap);
            caxis([vmin vmax]);
            cb = colorbar; ylabel(cb,'Correlation');

            xticks(1:n); yticks(1:n);
            xticklabels(labels); yticklabels(labels);
            xtickangle(45);
            for r=1:n, for c=1:n
                if ~isnan(M(r,c))
                    text(c,r,sprintf('%.3f',M(r,c)), ...
                        'HorizontalAlignment','center','FontSize',9,'Color','k');
                end
            end, end
            box on;
            title('Average Per-File Correlation (US 60Hz Experiment)');
            xlabel('Trace Type'); ylabel('Trace Type');

            saveas(fig, fullfile(results_dir,'avg_file_heatmap_US.svg'));
            fprintf('  Saved avg_file_heatmap_US.svg\n');
            fprintf(fileID,'  Saved avg_file_heatmap_US.svg\n');
            close(fig);
        else
            fprintf('  No US data found, skipping.\n');
            fprintf(fileID,'  No US data found, skipping.\n');
        end
    catch ME
        fprintf('  Error generating US average per-file heatmap: %s\n', ME.message);
        fprintf(fileID,'  Error generating US average per-file heatmap: %s\n', ME.message);
    end

    fprintf('STEP 2A Complete\n');
    fprintf(fileID,'STEP 2A Complete\n');
end

function generate_averaged_grid_freq_heatmap(all_results, results_dir, fileID)

    fprintf('\n\nSTEP 2B: Summarized Correlation Matrix Heatmap\n');
    fprintf(fileID, '\n\nSTEP 2B: Summarized Correlation Matrix Heatmap\n');

    vmin  = -0.06;
    vmax  =  1.0;
    gamma = 0.35;
    r_center = (vmin + vmax)/2;
    t_center = (r_center - vmin) / (vmax - vmin);
    t_center_g = t_center ^ gamma;
    numColors = 2048;
    cmap = setup_diverging_colormap(0, 1, t_center_g, numColors);

    all_grid = unique([{all_results.Grid_A}, {all_results.Grid_B}]);

    %% DE grids
    try
        de_names = all_grid(~cellfun('isempty', strfind(all_grid, 'DE_')));
        if ~isempty(de_names)
            n = numel(de_names);
            labels = cell(1,n);
            for i = 1:n
                labels{i} = get_corrmat_label_grid(de_names{i});
            end

            M = nan(n);
            for i = 1:n
                M(i,i) = 1.0;
                for j = 1:(i-1)
                    gi = de_names{i};
                    gj = de_names{j};
                    idx1 = find(strcmp({all_results.Grid_A},gi) & strcmp({all_results.Grid_B},gj));
                    idx2 = find(strcmp({all_results.Grid_A},gj) & strcmp({all_results.Grid_B},gi));
                    idx = [idx1 idx2];
                    if ~isempty(idx)
                        M(i,j) = mean([all_results(idx).Correlation]);
                    end
                end
            end

            T  = (M - vmin) / (vmax - vmin);
            T  = max(min(T,1),0);
            Tg = T .^ gamma;

            fig = figure('Name','DE Summary Heatmap', ...
                         'Visible','off','Position',[100 100 800 650]);
            imagesc(Tg, 'AlphaData', ~isnan(M));
            set(gca,'Color',[1 1 1]);
            axis equal tight;
            colormap(cmap);
            caxis([0 1]);

            cb = colorbar;
            nt = 6;
            ticks = linspace(0,1,nt);
            cb.Ticks = ticks;
            true_vals = vmin + (ticks .^ (1/gamma)) * (vmax - vmin);
            cb.TickLabels = arrayfun(@(x) sprintf('%.2f',x), true_vals, 'UniformOutput',false);
            ylabel(cb,'Correlation');

            xticks(1:n); yticks(1:n);
            xticklabels(labels); yticklabels(labels);
            xtickangle(45);

            for r = 1:n
                for c = 1:n
                    if ~isnan(M(r,c))
                        text(c, r, sprintf('%.3f', M(r,c)), ...
                             'HorizontalAlignment','center', ...
                             'FontSize',9,'Color','k');
                    end
                end
            end

            box on;
            title('Summary Correlation Matrix (DE 50Hz Experiment)');
            xlabel('Grid'); ylabel('Grid');

            fname = fullfile(results_dir, 'summary_correlation_heatmap_DE.svg');
            saveas(fig, fname);
            fprintf('  Saved %s\n', fname);
            fprintf(fileID,'  Saved %s\n', fname);
            close(fig);
        else
            fprintf('  No DE grid data found, skipping DE heatmap.\n');
            fprintf(fileID,'  No DE grid data found, skipping DE heatmap.\n');
        end
    catch ME
        fprintf('  Error generating DE summary heatmap: %s\n', ME.message);
        fprintf(fileID,'  Error generating DE summary heatmap: %s\n', ME.message);
    end

    %% US grids
    try
        us_names = all_grid(~cellfun('isempty', strfind(all_grid, 'US_')));
        if ~isempty(us_names)
            n = numel(us_names);
            labels = cell(1,n);
            for i = 1:n
                labels{i} = get_corrmat_label_grid(us_names{i});
            end

            M = nan(n);
            for i = 1:n
                M(i,i) = 1.0;
                for j = 1:(i-1)
                    gi = us_names{i};
                    gj = us_names{j};
                    idx1 = find(strcmp({all_results.Grid_A},gi) & strcmp({all_results.Grid_B},gj));
                    idx2 = find(strcmp({all_results.Grid_A},gj) & strcmp({all_results.Grid_B},gi));
                    idx = [idx1 idx2];
                    if ~isempty(idx)
                        M(i,j) = mean([all_results(idx).Correlation]);
                    end
                end
            end

            T  = (M - vmin) / (vmax - vmin);
            T  = max(min(T,1),0);
            Tg = T .^ gamma;

            fig = figure('Name','US Summary Heatmap', ...
                         'Visible','off','Position',[100 100 800 650]);
            imagesc(Tg, 'AlphaData', ~isnan(M));
            set(gca,'Color',[1 1 1]);
            axis equal tight;
            colormap(cmap);
            caxis([0 1]);

            cb = colorbar;
            nt = 6;
            ticks = linspace(0,1,nt);
            cb.Ticks = ticks;
            true_vals = vmin + (ticks .^ (1/gamma)) * (vmax - vmin);
            cb.TickLabels = arrayfun(@(x) sprintf('%.2f',x), true_vals, 'UniformOutput',false);
            ylabel(cb,'Correlation');

            xticks(1:n); yticks(1:n);
            xticklabels(labels); yticklabels(labels);
            xtickangle(45);

            for r = 1:n
                for c = 1:n
                    if ~isnan(M(r,c))
                        text(c, r, sprintf('%.3f', M(r,c)), ...
                             'HorizontalAlignment','center', ...
                             'FontSize',9,'Color','k');
                    end
                end
            end

            box on;
            title('Summary Correlation Matrix (US 60Hz Experiment)');
            xlabel('Grid'); ylabel('Grid');

            fname = fullfile(results_dir, 'summary_correlation_heatmap_US.svg');
            saveas(fig, fname);
            fprintf('  Saved %s\n', fname);
            fprintf(fileID,'  Saved %s\n', fname);
            close(fig);
        else
            fprintf('  No US grid data found, skipping US heatmap.\n');
            fprintf(fileID,'  No US grid data found, skipping US heatmap.\n');
        end
    catch ME
        fprintf('  Error generating US summary heatmap: %s\n', ME.message);
        fprintf(fileID,'  Error generating US summary heatmap: %s\n', ME.message);
    end

    fprintf('STEP 2B Complete\n');
    fprintf(fileID, 'STEP 2B Complete\n');
end

function cmap = setup_diverging_colormap(vmin, vmax, vcenter, numColors)
    % Diverging colormap with extra resolution near the upper end.
    % - vmin maps to red
    % - vcenter maps to white
    % - vmax maps to deep blue
    % - Majority of colors concentrated between 0.9 and 1.0 so that
    %   values like 0.998 vs 1.000 get visibly different colors.

    if nargin < 4
        numColors = 2048;  % high resolution
    end
    if nargin < 3 || isempty(vcenter)
        vcenter = (vmin + vmax) / 2;
    end

    % Anchor colors
    neg_dark = [229,115,115]/255;   % red
    neu_color = [1,1,1];            % white
    pos_dark = [21,101,192]/255;    % deep blue

    span = vmax - vmin;

    % Piecewise value mapping:
    %  - First 20% of colormap:  vmin -> 0.9   (compressed)
    %  - Remaining 80%:          0.9  -> vmax (expanded; high resolution)
    split_val = 0.9;
    split_alpha = 0.20;  % fraction of colormap for [vmin, split_val]

    idx = linspace(0,1,numColors);  % normalized colormap index in [0,1]
    val = zeros(size(idx));

    low_mask = idx <= split_alpha;
    hi_mask  = ~low_mask;

    if any(low_mask)
        val(low_mask) = vmin + (split_val - vmin) * (idx(low_mask) / split_alpha);
    end
    if any(hi_mask)
        val(hi_mask) = split_val + (vmax - split_val) * ((idx(hi_mask) - split_alpha) / (1 - split_alpha));
    end

    % Now map each "val" to RGB via diverging scheme
    cmap = zeros(numColors, 3);
    for k = 1:numColors
        if val(k) <= vcenter
            t = (val(k) - vmin) / max(vcenter - vmin, eps);
            t = max(min(t,1),0);
            cmap(k,:) = neg_dark + t * (neu_color - neg_dark);
        else
            t = (val(k) - vcenter) / max(vmax - vcenter, eps);
            t = max(min(t,1),0);
            cmap(k,:) = neu_color + t * (pos_dark - neu_color);
        end
    end
end

function enable_t_test_analysis(intra_grid_correlations_z, inter_grid_correlations_z, fileID)
    fprintf('\n\nSTEP 3: Formal Hypothesis Testing\n');
    fprintf(fileID, '\n\nSTEP 3: Formal Hypothesis Testing\n');
    
    try
        if isempty(intra_grid_correlations_z) || isempty(inter_grid_correlations_z)
            fprintf('  Skipping t-test: One or both correlation groups are empty.\n');
            fprintf(fileID, '  Skipping t-test: One or both correlation groups are empty.\n');
            return;
        end
        
        fprintf('  Performing two-sample t-test (intra-grid > inter-grid)...\n');
        fprintf(fileID, '  Performing two-sample t-test (intra-grid > inter-grid)...\n');

        z_intra = intra_grid_correlations_z;      
        z_inter = inter_grid_correlations_z;
    
        [~, p_value, ~, stats] = ttest2(z_intra, z_inter, 'Vartype', 'unequal', 'Tail', 'right');
    
        p_string     = sprintf('p = %.4e', p_value);
        N_intra      = length(z_intra);
        N_inter      = length(z_inter);
        z_intra_mean = mean(z_intra);
        z_inter_mean = mean(z_inter);
        z_intra_std  = std(z_intra);
        z_inter_std  = std(z_inter);
        
        std_pooled = sqrt( ...
            ((N_intra-1)*z_intra_std^2 + (N_inter-1)*z_inter_std^2) / ...
            (N_intra + N_inter - 2) );
        cohens_d = (z_intra_mean - z_inter_mean) / std_pooled;
    
        fprintf('  t-test on Fisher-transformed z-values:\n');
        fprintf(fileID, '  t-test on Fisher-transformed z-values:\n');
        fprintf('    t-statistic: %.4f\n', stats.tstat);
        fprintf(fileID, '    t-statistic: %.4f\n', stats.tstat);
        fprintf('    Degrees of Freedom (df): %.4f\n', stats.df);
        fprintf(fileID, '    Degrees of Freedom (df): %.4f\n', stats.df);
        fprintf('    p-value: %s\n', p_string);
        fprintf(fileID, '    p-value: %s\n', p_string);    
        fprintf('    Effect Size (Cohen''s d): %.4f\n', cohens_d);
        fprintf(fileID, '    Effect Size (Cohen''s d): %.4f\n', cohens_d);
        fprintf('\n');
        fprintf(fileID, '\n');      
        
    catch ME_ttest
        fprintf('  Error during formal hypothesis testing: %s\n', ME_ttest.message);
        fprintf(fileID, '  Error during formal hypothesis testing: %s\n', ME_ttest.message);
    end
    
    fprintf('STEP 3 Complete\n');
    fprintf(fileID, 'STEP 3 Complete\n');
end

function plot_z_score_distributions(z_stats, results_dir, fileID)
% Plot parametric separation between intra- and inter-grid scores
% using normal fits on Fisher z'-transformed correlations.

    % ---------- Summary stats ----------
    mu_intra = z_stats.intra.mean;
    sd_intra = z_stats.intra.std;
    N_intra  = z_stats.intra.N;

    mu_inter = z_stats.inter.mean;
    sd_inter = z_stats.inter.std;
    N_inter  = z_stats.inter.N;

    se_intra = sd_intra / sqrt(N_intra);
    se_inter = sd_inter / sqrt(N_inter);
    ci_intra = [mu_intra - 1.96*se_intra, mu_intra + 1.96*se_intra];
    ci_inter = [mu_inter - 1.96*se_inter, mu_inter + 1.96*se_inter];

    % ---------- Range & PDFs ----------
    lo = min(mu_inter - 4*sd_inter, mu_intra - 4*sd_intra);
    hi = max(mu_inter + 4*sd_inter, mu_intra + 4*sd_intra);
    x  = linspace(lo, hi, 2000);

    pdf_intra = normpdf(x, mu_intra, sd_intra);
    pdf_inter = normpdf(x, mu_inter, sd_inter);

    % ---------- Overlap ----------
    overlap_curve = min(pdf_intra, pdf_inter);
    overlap_area  = trapz(x, overlap_curve);
    overlap_pct   = overlap_area * 100;

    fprintf('    > Overlap Area between distributions: %.6f (%.4f%%)\n', ...
            overlap_area, overlap_pct);
    fprintf(fileID, ...
            '    > Overlap Area between distributions: %.6f (%.4f%%)\n', ...
            overlap_area, overlap_pct);

    % ---------- Parametric EER threshold ----------
    tau_star = NaN;
    eer_val  = NaN;
    if all(isfinite([mu_inter, sd_inter, mu_intra, sd_intra])) && ...
       sd_inter > 0 && sd_intra > 0

        tau_star = (sd_inter*mu_intra + sd_intra*mu_inter) / (sd_inter + sd_intra);

        z0 = (tau_star - mu_inter) / sd_inter;
        z1 = (tau_star - mu_intra) / sd_intra;
        fpr_eer = 1 - normcdf(z0);
        fnr_eer = normcdf(z1);
        eer_val = 0.5 * (fpr_eer + fnr_eer);
    end

    % ---------- Colors ----------
    col_intra = [0.0000 0.4470 0.7410];  % blue
    col_inter = [0.8500 0.3250 0.0980];  % orange
    col_ovlp_main  = [0.85 0.85 0.95];   % light for main
    col_ovlp_inset = [0.75 0.75 0.75];   % gray for inset

    % ---------- Figure & main axes ----------
    fig = figure('Name', 'Z-Transformed Distributions', ...
                 'Visible', 'off', ...
                 'Position', [100, 100, 1000, 600]);

    % Shrink main axes width to leave space on the right for inset
    ax = axes('Parent', fig, 'Position', [0.08 0.18 0.60 0.75]);
    hold(ax, 'on');

    % PDFs (primary)
    h_intra = plot(ax, x, pdf_intra, 'Color', col_intra, 'LineWidth', 2.2);
    h_inter = plot(ax, x, pdf_inter, 'Color', col_inter, 'LineWidth', 2.2);

    % Overlap (light)
    h_ovlp = fill(ax, x, overlap_curve, col_ovlp_main, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.7);

    % Mean lines
    yl = ylim(ax);
    h_mean_intra = line(ax, [mu_intra mu_intra], yl, ...
                        'Color', col_intra, 'LineStyle', '--', 'LineWidth', 1.2);
    h_mean_inter = line(ax, [mu_inter mu_inter], yl, ...
                        'Color', col_inter, 'LineStyle', '--', 'LineWidth', 1.2);

    % EER threshold (if defined)
    if ~isnan(tau_star) && ~isnan(eer_val)
        h_tau = line(ax, [tau_star tau_star], yl, ...
                     'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 1.5);
    else
        h_tau = gobjects(0);
    end
    ylim(ax, yl);

    % ---------- Axes styling & ticks ----------
    xlim(ax, [lo hi]);
    grid(ax, 'on');
    ax.GridAlpha = 0.15;
    ax.Box = 'on';

    % More ticks:
    x_step = 0.5;
    xticks(ax, ceil(lo/x_step)*x_step : x_step : floor(hi/x_step)*x_step);

    max_pdf = max([pdf_intra, pdf_inter]);
    y_step = 0.1;
    yticks(ax, 0 : y_step : ceil(max_pdf / y_step) * y_step);

    xlabel(ax, 'Z-Score (z'')');
    ylabel(ax, 'Probability Density');
    title(ax, 'Separation of Intra-Grid vs Inter-Grid Scores (Fisher z'')');

    set(ax, 'FontName', 'Helvetica', 'FontSize', 9);

    % ---------- Legend at bottom ----------
    leg_entries = {
        sprintf('Intra (\\mu=%.2f, \\sigma=%.2f, 95%%%% CI [%.2f, %.2f])', ...
                mu_intra, sd_intra, ci_intra(1), ci_intra(2)), ...
        sprintf('Inter (\\mu=%.2f, \\sigma=%.2f, 95%%%% CI [%.2f, %.2f])', ...
                mu_inter, sd_inter, ci_inter(1), ci_inter(2)), ...
        sprintf('Overlap = %.4f%%%%', overlap_pct)
    };

    h_for_legend = [h_intra, h_inter, h_ovlp];

    if ~isempty(h_tau)
        leg_entries{end+1} = sprintf('\\tau^* (EER \\approx %.2e)', eer_val);
        h_for_legend(end+1) = h_tau;
    end

    lgd = legend(ax, h_for_legend, leg_entries, ...
                 'Location', 'southoutside', ...
                 'Orientation', 'vertical', ...
                 'Box', 'on', ...
                 'FontSize', 8);
    % Slightly tighten legend box
    lgd.Position(2) = 0.02;

        % ---------- Inset: zoomed overlap on the right ----------
    try
        ax_inset = axes('Parent', fig, 'Position', [0.72 0.45 0.22 0.40]);
        hold(ax_inset, 'on');

        % PDFs in inset
        plot(ax_inset, x, pdf_intra, 'Color', col_intra, 'LineWidth', 1.0);
        plot(ax_inset, x, pdf_inter, 'Color', col_inter, 'LineWidth', 1.0);

        % Gray overlap shading
        fill(ax_inset, x, overlap_curve, col_ovlp_inset, ...
             'EdgeColor', 'none', 'FaceAlpha', 0.7);

        % Center zoom around max overlap
        [~, idx_max_ovlp] = max(overlap_curve);
        x_c = x(idx_max_ovlp);
        w   = max(sd_intra, sd_inter) * 3;   % span around crossing

        x_min_z = max(lo, x_c - w/2);
        x_max_z = min(hi, x_c + w/2);
        if x_max_z <= x_min_z
            x_min_z = x_c - 0.5;
            x_max_z = x_c + 0.5;
        end
        xlim(ax_inset, [x_min_z, x_max_z]);

        % Y-limits based on local overlap
        local_idx   = (x >= x_min_z) & (x <= x_max_z);
        max_pdf_z   = max(overlap_curve(local_idx));
        if isempty(max_pdf_z) || max_pdf_z <= 0
            max_pdf_z = max(pdf_intra(local_idx) + pdf_inter(local_idx));
        end
        if isempty(max_pdf_z) || max_pdf_z <= 0
            max_pdf_z = max_pdf / 5;
        end
        ylim(ax_inset, [0, max_pdf_z * 1.2]);

        % High-resolution x-ticks (? 6–10 ticks, step snapped to 0.01)
        desired_ticks = 8;
        raw_step = (x_max_z - x_min_z) / max(desired_ticks - 1, 1);
        step = max(0.01, round(raw_step / 0.01) * 0.01);  % 0.01, 0.02, 0.03, ...
        xt = x_min_z : step : x_max_z;
        if numel(xt) < 4
            step = step / 2;
            xt = x_min_z : step : x_max_z;
        end
        xt = round(xt, 3);  % avoid ugly floating labels

        xticks(ax_inset, xt);
        xticklabels(ax_inset, arrayfun(@(v) sprintf('%.2f', v), xt, 'UniformOutput', false));

        title(ax_inset, 'Overlap (Zoomed)', 'FontSize', 8);
        xlabel(ax_inset, 'z''', 'FontSize', 7);
        ylabel(ax_inset, 'Density', 'FontSize', 7);

        set(ax_inset, 'FontSize', 7, 'Box', 'on');
        grid(ax_inset, 'on');
        ax_inset.GridAlpha = 0.15;

    catch ME_zoom
        fprintf('    Warning: Could not generate zoomed overlap plot: %s\n', ME_zoom.message);
        fprintf(fileID, '    Warning: Could not generate zoomed overlap plot: %s\n', ME_zoom.message);
    end


    % ---------- Save as SVG ----------
    out_name = fullfile(results_dir, 'z_score_distributions.svg');
    try
        saveas(fig, out_name);
        fprintf('    Saved %s\n', out_name);
        fprintf(fileID, '    Saved %s\n', out_name);
    catch ME
        fprintf('    Error saving SVG: %s\n', ME.message);
        fprintf(fileID, '    Error saving SVG: %s\n', ME.message);
    end

    close(fig);
end

function identity = get_file_identity(filename, baseDir, nominal_freq_arr)
    identity.filename = filename;
    
    if ~isempty(strfind(baseDir, 'DE_50'))
        if ~isempty(strfind(filename, 'dresden'))
            identity.grid = 'DE_50Hz';
            identity.freq = 50;
        else
            identity.grid = 'DE_Control_60Hz'; 
            identity.freq = 60;
        end
    elseif ~isempty(strfind(baseDir, 'US_60'))
        identity.freq = 60;
        if ~isempty(strfind(filename, 'egrid'))
            identity.grid = 'US_EGrid_60Hz';
        elseif ~isempty(strfind(filename, 'tgrid'))
            identity.grid = 'US_TGrid_60Hz';
        elseif ~isempty(strfind(filename, 'wgrid'))
            identity.grid = 'US_WGrid_60Hz';
        else
            identity.grid = 'US_Unknown_60Hz';
        end
    else
        identity.grid = 'US_EGrid_60Hz';
        identity.freq = 60; 
    end
    
    if identity.freq == 50
        params.nominal_freq = nominal_freq_arr(1);
        params.plot_title   = sprintf('%s (50 Hz)', filename);
    else
        params.nominal_freq = nominal_freq_arr(2);
        params.plot_title   = sprintf('%s (60 Hz)', filename);
    end
    
    params.harmonics_arr           = (1:7) * params.nominal_freq;
    params.est_freq                = params.harmonics_arr(1);
    params.est_method              = 1;
    params.est_spec_comb_harmonics = params.harmonics_arr;
    identity.params = params;
end

function label = get_corrmat_label(filename)
    if ~isempty(strfind(filename, 'fpga_em_trace_dc_egrid_citya_lab'))
        label = 'EM'; 
    elseif ~isempty(strfind(filename, 'mains_pow_trace_ac_egrid_citya_lab'))
        label = 'AR1'; 
    elseif ~isempty(strfind(filename, 'mains_pow_trace_ac_egrid_citya_home'))
        label = 'AR2';
    elseif ~isempty(strfind(filename, 'mains_pow_trace_ac_egrid_worcester'))
        label = 'AR3'; 
    elseif ~isempty(strfind(filename, 'mains_pow_trace_ac_tgrid_richardson'))
        label = 'AR4'; 
    elseif ~isempty(strfind(filename, 'mains_pow_trace_ac_wgrid_tucson'))
        label = 'AR5'; 
    elseif ~isempty(strfind(filename, 'fpga_em_trace_dc_citya_lab'))
        label = 'EM_DE'; 
    elseif ~isempty(strfind(filename, 'mains_pow_trace_ac_citya_lab'))
        label = 'AR_DE'; 
    elseif ~isempty(strfind(filename, 'mains_pow_trace_ac_dresden'))
        label = 'AR_Dres'; 
    else
        label = strrep(filename, '_', '-'); 
    end
end

function label = get_corrmat_label_grid(grid_name)
    if ~isempty(strfind(grid_name, 'DE_50Hz'))
        label = 'Germany 50Hz';
    elseif ~isempty(strfind(grid_name, 'DE_Control_60Hz'))
        label = 'US East 60Hz';
    elseif ~isempty(strfind(grid_name, 'US_EGrid_60Hz'))
        label = 'US East 60Hz';
    elseif ~isempty(strfind(grid_name, 'US_TGrid_60Hz'))
        label = 'US Texas 60Hz';
    elseif ~isempty(strfind(grid_name, 'US_WGrid_60Hz'))
        label = 'US West 60Hz';
    elseif ~isempty(strfind(grid_name, 'US_Unknown_60Hz'))
        label = 'Unknown 60Hz';
    else
        label = strrep(grid_name, '_', '-');
    end
end

function out = iif(condition, trueValue, falseValue)
    if condition
        out = trueValue;
    else
        out = falseValue;
    end
end

function add_heatmap_caption(fig_handle, caption_str)
    try
        fig_pos = get(fig_handle, 'Position');
        original_height = fig_pos(4);
        caption_height  = 40;
        
        fig_pos(4) = original_height + caption_height;
        fig_pos(2) = max(1, fig_pos(2) - caption_height);
        set(fig_handle, 'Position', fig_pos);
        
        ax = findobj(fig_handle, 'Type', 'Axes');
        if ~isempty(ax)
            main_ax = [];
            for a_idx = 1:length(ax)
                if ~isa(ax(a_idx), 'matlab.graphics.illustration.ColorBar')
                    main_ax = ax(a_idx);
                    break;
                end
            end
            
            if ~isempty(main_ax)
                ax_pos = get(main_ax, 'Position');
                bottom_margin = caption_height / fig_pos(4);
                ax_pos(2) = ax_pos(2) + bottom_margin;
                ax_pos(4) = ax_pos(4) - bottom_margin;
                set(main_ax, 'Position', ax_pos);
            end
        end
        
        annotation(fig_handle, 'textbox', [0.1, 0.01, 0.8, 0.08], 'String', caption_str, ...
          'EdgeColor','none','HorizontalAlignment','left','Interpreter','none', ...
          'FontSize', 8, 'FitBoxToText', 'off');
    catch
    end
end

function caption_str = generate_heatmap_caption(labels)
    definitions = {};
    if any(strcmp(labels, 'EM'))
        definitions{end+1} = 'EM = EM Sens (EG, CtyA Lab)';
    end
    if any(strcmp(labels, 'AR1'))
        definitions{end+1} = 'AR1 = AC Ref (EG, CtyA Lab)';
    end
    if any(strcmp(labels, 'AR2'))
        definitions{end+1} = 'AR2 = AC Ref (EG, CtyA Home)';
    end
    if any(strcmp(labels, 'AR3'))
        definitions{end+1} = 'AR3 = AC Ref (EG, Worc.)';
    end
    if any(strcmp(labels, 'AR4'))
        definitions{end+1} = 'AR4 = AC Ref (TG, Rich.)';
    end
    if any(strcmp(labels, 'AR5'))
        definitions{end+1} = 'AR5 = AC Ref (WG, Tucs.)';
    end
    if any(strcmp(labels, 'EM_DE'))
        definitions{end+1} = 'EM_DE = EM Sens (CtyA Lab, 60Hz)';
    end
    if any(strcmp(labels, 'AR_DE'))
        definitions{end+1} = 'AR_DE = AC Ref (CtyA Lab, 60Hz)';
    end
    if any(strcmp(labels, 'AR_Dres'))
        definitions{end+1} = 'AR_Dres = AC Ref (Dres., 50Hz)';
    end
    
    caption_str = strjoin(definitions, '; ');
end

function eer_stats = compute_parametric_eer_gaussian(z_stats, fileID)
    % Classes:
    %   inter ~ N(mu0, s0)      (impostor / different-grid)
    %   intra ~ N(mu1, s1)      (genuine / same-grid)
    % Decision rule (on Fisher z'): accept "at expected site" if z >= tau.
    %
    % EER (unequal-variance) threshold:
    %   tau* solves FPR(tau) = FNR(tau) and equals:
    %   tau* = (mu0*s1 + mu1*s0) / (s0 + s1)
    %   EER  = 1 - Phi((tau* - mu0)/s0) = Phi((tau* - mu1)/s1)
    %
    % Fixed-FPR alpha:
    %   tau(alpha) = mu0 + s0 * Phi^{-1}(1 - alpha)
    %   FNR(alpha) = Phi((tau(alpha) - mu1)/s1)
    %
    % Fixed-FNR beta:
    %   tau(beta)  = mu1 + s1 * Phi^{-1}(beta)
    %   FPR(beta)  = 1 - Phi((tau(beta) - mu0)/s0)

    % --- Means & stds from data
    mu0 = z_stats.inter.mean;
    s0  = z_stats.inter.std;
    mu1 = z_stats.intra.mean;
    s1  = z_stats.intra.std;

    % --- Adjusted std for temporal variation (apply to intra by default)
    % Prefer a field provided by upstream code; else use your known value.
    if isfield(z_stats.intra, 'std_adj') && ~isempty(z_stats.intra.std_adj)
        s1_adj = z_stats.intra.std_adj;
    else
        s1_adj = 0.6366;  % <-- from temporal-variation std for intra (same-grid)
    end

    % Sanity check
    if any(isnan([mu0,mu1,s0,s1,s1_adj])) || any([s0,s1,s1_adj] <= 0)
        tau_star    = NaN; eer = NaN; fpr_eer = NaN; fnr_eer = NaN;
        tau_star_ad = NaN; eer_adj = NaN;
    else
        % --------- ORIGINAL (unadjusted) EER ----------
        tau_star = (mu0*s1 + mu1*s0) / (s0 + s1);
        z0 = (tau_star - mu0)/s0;
        z1 = (tau_star - mu1)/s1;
        fpr_eer = 1 - normcdf(z0);
        fnr_eer = normcdf(z1);
        %eer     = 0.5*(fpr_eer + fnr_eer);  % numerically identical here
        eer = normcdf((mu0 - mu1)/(s0 + s1));

        % --------- ADJUSTED EER (replace intra std by s1_adj) ----------
        tau_star_ad = (mu0*s1_adj + mu1*s0) / (s0 + s1_adj);
        z0_ad = (tau_star_ad - mu0)/s0;
        z1_ad = (tau_star_ad - mu1)/s1_adj;
        fpr_eer_ad = 1 - normcdf(z0_ad);
        fnr_eer_ad = normcdf(z1_ad);
        %eer_adj    = 0.5*(fpr_eer_ad + fnr_eer_ad);
        eer_adj = normcdf((mu0 - mu1)/(s0 + s1_adj));
    end

    % --------- Operating points for fixed FPR or FNR ----------
    % Targets you care about (edit as needed)
    alpha_targets = [1e-6];  % fixed FPR = alpha
    beta_targets  = [1e-6];  % fixed FNR = beta

    % ORIGINAL (uses s0, s1)
    ops_fixed_fpr = arrayfun(@(a) fixed_fpr_point(a, mu0, s0, mu1, s1), alpha_targets);
    ops_fixed_fnr = arrayfun(@(b) fixed_fnr_point(b, mu0, s0, mu1, s1), beta_targets);

    % ADJUSTED (uses s0, s1_adj)
    ops_fixed_fpr_adj = arrayfun(@(a) fixed_fpr_point(a, mu0, s0, mu1, s1_adj), alpha_targets);
    ops_fixed_fnr_adj = arrayfun(@(b) fixed_fnr_point(b, mu0, s0, mu1, s1_adj), beta_targets);

    % --------- Package outputs ----------
    eer_stats = struct();
    eer_stats.mu0 = mu0; eer_stats.sigma0 = s0;
    eer_stats.mu1 = mu1; eer_stats.sigma1 = s1;
    eer_stats.sigma1_adj = s1_adj;

    eer_stats.EER = struct('tau', tau_star, 'FPR', fpr_eer, 'FNR', fnr_eer, 'rate', eer);
    eer_stats.EER_ADJ = struct('tau', tau_star_ad, 'FPR', fpr_eer_ad, 'FNR', fnr_eer_ad, 'rate', eer_adj);

    % Fixed-FPR tables
    eer_stats.fixed_FPR = struct( ...
        'targets', alpha_targets, ...
        'original', ops_fixed_fpr, ...
        'adjusted', ops_fixed_fpr_adj);

    % Fixed-FNR tables
    eer_stats.fixed_FNR = struct( ...
        'targets', beta_targets, ...
        'original', ops_fixed_fnr, ...
        'adjusted', ops_fixed_fnr_adj);

    % --------- Logging ----------
    fprintf('\nSTEP 4: Parametric EER (Gaussian on Fisher z'')\n');
    fprintf('  mu0 (inter) = %.4f, sigma0 = %.4f\n', mu0, s0);
    fprintf('  mu1 (intra) = %.4f, sigma1 = %.4f\n', mu1, s1);
    if ~isnan(eer_stats.EER.rate)
        fprintf('  tau* (orig) = %.4f | EER ~= %.4e  [FPR=%.3e, FNR=%.3e]\n', ...
            eer_stats.EER.tau, eer_stats.EER.rate, eer_stats.EER.FPR, eer_stats.EER.FNR);
        fprintf('  tau* (ADJ)  = %.4f | EER_ADJ ~= %.4e  [FPR=%.3e, FNR=%.3e]  (sigma1_adj=%.4f)\n', ...
            eer_stats.EER_ADJ.tau, eer_stats.EER_ADJ.rate, ...
            eer_stats.EER_ADJ.FPR, eer_stats.EER_ADJ.FNR, s1_adj);

        for k = 1:numel(alpha_targets)
            a = alpha_targets(k);
            op  = ops_fixed_fpr(k);
            opA = ops_fixed_fpr_adj(k);
            fprintf('  Fix FPR=%.1e:  tau=%.4f => FNR(orig)=%.3e | tau_ADJ=%.4f => FNR(ADJ)=%.3e\n', ...
                a, op.tau, op.FNR, opA.tau, opA.FNR);
        end
        for k = 1:numel(beta_targets)
            b = beta_targets(k);
            op  = ops_fixed_fnr(k);
            opA = ops_fixed_fnr_adj(k);
            fprintf('  Fix FNR=%.1e:  tau=%.4f => FPR(orig)=%.3e | tau_ADJ=%.4f => FPR(ADJ)=%.3e\n', ...
                b, op.tau, op.FPR, opA.tau, opA.FPR);
        end
    else
        fprintf('  EER = NaN (insufficient or degenerate data)\n');
    end

    if fileID > 1
        fprintf(fileID, '\nSTEP 4: Parametric EER (Gaussian on Fisher z'')\n');
        fprintf(fileID, '  mu0 (inter) = %.4f, sigma0 = %.4f\n', mu0, s0);
        fprintf(fileID, '  mu1 (intra) = %.4f, sigma1 = %.4f (sigma1_adj=%.4f)\n', mu1, s1, s1_adj);
        if ~isnan(eer_stats.EER.rate)
            fprintf(fileID, '  tau* (orig)=%.4f | EER ~= %.4e [FPR=%.3e, FNR=%.3e]\n', ...
                eer_stats.EER.tau, eer_stats.EER.rate, eer_stats.EER.FPR, eer_stats.EER.FNR);
            fprintf(fileID, '  tau* (ADJ)=%.4f | EER_ADJ ~= %.4e [FPR=%.3e, FNR=%.3e]\n', ...
                eer_stats.EER_ADJ.tau, eer_stats.EER_ADJ.rate, ...
                eer_stats.EER_ADJ.FPR, eer_stats.EER_ADJ.FNR);

            for k = 1:numel(alpha_targets)
                a = alpha_targets(k);
                op  = ops_fixed_fpr(k);
                opA = ops_fixed_fpr_adj(k);
                fprintf(fileID, '  Fix FPR=%.1e: tau=%.4f => FNR(orig)=%.3e | tau_ADJ=%.4f => FNR(ADJ)=%.3e\n', ...
                    a, op.tau, op.FNR, opA.tau, opA.FNR);
            end
            for k = 1:numel(beta_targets)
                b = beta_targets(k);
                op  = ops_fixed_fnr(k);
                opA = ops_fixed_fnr_adj(k);
                fprintf(fileID, '  Fix FNR=%.1e: tau=%.4f => FPR(orig)=%.3e | tau_ADJ=%.4f => FPR(ADJ)=%.3e\n', ...
                    b, op.tau, op.FPR, opA.tau, opA.FPR);
            end
        else
            fprintf(fileID, '  EER = NaN (insufficient or degenerate data)\n');
        end
    end
end

function out = fixed_fpr_point(alpha, mu0, s0, mu1, s1)
    % tau that yields FPR = alpha
    tau = mu0 + s0 * norminv(1 - alpha);
    fnr = normcdf((tau - mu1)/s1);
    out = struct('tau', tau, 'FPR', alpha, 'FNR', fnr);
end

function out = fixed_fnr_point(beta, mu0, s0, mu1, s1)
    % tau that yields FNR = beta
    tau = mu1 + s1 * norminv(beta);
    fpr = 1 - normcdf((tau - mu0)/s0);
    out = struct('tau', tau, 'FNR', beta, 'FPR', fpr);
end

