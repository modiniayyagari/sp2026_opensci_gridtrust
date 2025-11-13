%% ENF Server-Room Robustness (Experiment 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   Evaluates robustness of ambient-EM ENF extraction in a high-density
%   server-room environment. Tests whether sensor position on the server PSU
%   and day-of-week conditions materially affect the correlation between the
%   sensed EM ENF and the mains reference. Produces summary statistics and
%   model-based mean ±95% CI (from Fisher z') using 1-way repeated-measures
%   ANOVA for Site (SBOX vs SPSU) and Day (MON vs WED).
%
% ANALYSIS OVERVIEW
% 1) INPUTS
%    • Root study directory with subfolders labeled by:
%        Site ? {SBOX, SPSU}, Day ? {MON, WED}, Time ? {T01…T05}.
%    • For each condition, two time-aligned WAV files:
%        - mains_pow_trace_ac.wav     (ground-truth mains reference)
%        - fpga_em_trace_dc.wav       (ambient EM sensed near the server PSU)
%
% 2) ENF EXTRACTION & CORRELATION
%    • STFT around 60-Hz harmonics (configurable nfft/frame/overlap).
%    • Weighted-harmonic ENF estimator per trace.
%    • Pearson correlation r for each sensed/reference pair.
%    • Fisher transform z = atanh(r) for parametric inference.
%
% 3) STATISTICS
%    • Descriptives on r and z (mean, SD, SEM, 95% CI).
%    • One-way RM-ANOVA on z for:
%        - Day: MON vs WED (subjects = Time, averaged over Site)
%        - Site: SBOX vs SPSU (subjects = Time, averaged over Day)
%      Assumption checks (Anderson–Darling normality, Levene homogeneity).
%      Partial ?_p² reported; with 2 levels per factor, GG = raw p.
%    • Back-transform fixed-effect estimates and CIs to r for interpretability.
%
% 4) VISUALIZATION & REPORTING
%    • Dot-and-whisker plots of mean r ±95% CI for Day and Site.
%    • Boxplots and histograms for r and z; z-distribution diagnostic plot.
%    • Fixed-width ANOVA summary tables printed to the log.
%
% 5) OUTPUTS
%    • Logs ? exp_logs/
%    • Figures ? exp_results/enf_analysis_top_reliability_stat_server/ (or config.results_dir)
%      All figures saved as .svg.
%
% Usage:
%   Set config.baseDir to the SRV_L root and ensure proc_enf_analysis.p is on the
%   MATLAB path. Run this script to generate the full analysis and figures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% 1) Configuration
config = struct();

%Main Path Configuration
config.baseDir = '../exp_inputs/SRV_L/';

%File Name Configuration
config.file_1_name_base = 'mains_pow_trace_ac'; % reference
config.file_2_name_base = 'fpga_em_trace_dc';   % sensed

% Expected labels
config.allowed_sites = {'SBOX','SPSU'};
config.allowed_days  = {'MON','WED'};
config.allowed_times = {'T01','T02','T03','T04','T05'};

% Logging / plotting flags
config.SCR_CFG_DISP_RAW_PAIR_WISE_CORR_VALS = false;
config.SCR_CFG_DISP_DESCRIPTIVE_STATS_RES   = true;
config.SCR_CFG_LIST_SUBJECT_KEYS            = false;
config.SCR_CFG_QQ_PLOTS                     = true;

% ENF extraction parameters
config.nominal_freq_arr = [50 60];      % [50Hz, 60Hz]
frame_size_arr          = (1:12)*1000;
config.frame_size       = frame_size_arr(8);  % 8000 ms
nfft_arr                = 2.^(10:20);
config.nfft             = nfft_arr(6);        % 2^15
overlap_size_arr        = 0:0.1:0.9;
config.overlap_size     = overlap_size_arr(1)*config.frame_size;

% Trace-specific parameters (for proc_enf_analysis)
nominal_freq = config.nominal_freq_arr(2); % 60 Hz
harmonics = (1:7) * nominal_freq;

config.harmonics_arr_1                        = harmonics;
config.trace_1_freq_est_method                = 1;
config.trace_1_est_freq                       = nominal_freq;
config.trace_1_freq_est_spec_comb_harmonics   = harmonics;
config.trace_1_plot_title                     = config.file_1_name_base;

config.harmonics_arr_2                        = harmonics;
config.trace_2_freq_est_method                = 1;
config.trace_2_est_freq                       = nominal_freq;
config.trace_2_freq_est_spec_comb_harmonics   = harmonics;
config.trace_2_plot_title                     = config.file_2_name_base;

%% 2) Setup logging/output folders
[fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders();
config.results_dir = results_dir;

fprintf('Starting ENF Server-Room Robustness Analysis. Log: %s\n', log_filename);
fprintf(fileID, 'Starting ENF Server-Room Robustness Analysis. Log: %s\n\n', log_filename);

%% 3) Scan folders and compute correlations
tbl = run_serverroom_scan(config, fileID);  % Correlation, FilePath, Site, Day, Time

if isempty(tbl)
    fprintf('No data found. Exiting.\n');
    fprintf(fileID, 'No data found. Exiting.\n');
    return;
end

all_r = tbl.Correlation;

%% 4) Descriptive stats + summary plots
[r_stats, z_data, z_stats] = calculate_basic_stats(all_r);

if config.SCR_CFG_DISP_DESCRIPTIVE_STATS_RES
    print_descriptive_stats_results(r_stats, z_stats, fileID);
    try
        plot_box_simple(all_r, r_stats, results_dir, fileID, 'r');
        plot_hist_simple(all_r, r_stats, results_dir, fileID, 'r');
        plot_box_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_hist_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_z_score_distributions(z_data, z_stats, results_dir, fileID);
    catch MEp
        fprintf('Plot error: %s\n', MEp.message);
        fprintf(fileID, 'Plot error: %s\n', MEp.message);
    end
else
    fprintf('\nSTEP 1 (Descriptive) skipped by config.\n');
    fprintf(fileID, '\nSTEP 1 (Descriptive) skipped by config.\n');
end

if config.SCR_CFG_DISP_RAW_PAIR_WISE_CORR_VALS
    print_raw_correlation_data_noID(tbl, all_r, fileID);
end

%% 5) [DAY] 1-way RM-ANOVA (MON vs WED), subjects = Time (averaged over Site)
fprintf('\n[DAY] 1-way RM-ANOVA (MON vs WED), subjects = Time\n');
fprintf(fileID, '\n[DAY] 1-way RM-ANOVA (MON vs WED), subjects = Time\n');
run_oneway_rm_day(tbl, config, fileID);

%% 6) [SITE] 1-way RM-ANOVA (SBOX vs SPSU), subjects = Time (averaged over Day)
fprintf('\n[SITE] 1-way RM-ANOVA (SBOX vs SPSU), subjects = Time\n');
fprintf(fileID, '\n[SITE] 1-way RM-ANOVA (SBOX vs SPSU), subjects = Time\n');
run_oneway_rm_site(tbl, config, fileID);

fprintf('\nDone. Results in: %s\n', results_dir);
fprintf(fileID, '\nDone. Results in: %s\n', results_dir);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedures for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders()
    log_dir     = 'exp_logs';
    results_dir = 'exp_results/enf_analysis_top_reliability_stat_server';
    
    if ~exist(log_dir,'dir'),     mkdir(log_dir);     end
    if ~exist(results_dir,'dir'), mkdir(results_dir); end
    
    [~, script_name] = fileparts(mfilename('fullpath'));
    if isempty(script_name), script_name = 'enf_analysis_top_reliability_stat_server'; end
    
    log_filename = fullfile(log_dir, [script_name '_log.txt']);
    if exist(log_filename,'file'), delete(log_filename); end
    
    fileID = fopen(log_filename,'w');
    if fileID == -1, error('Could not open %s', log_filename); end
    cleanupObj = onCleanup(@() fclose(fileID));
end

function tbl = run_serverroom_scan(config, fileID)
    % Build results table
    varTypes = {'double','string','string','string','string'};
    varNames = {'Correlation','FilePath','Site','Day','Time'};
    tbl = table('Size',[0,numel(varNames)], ...
                'VariableTypes',varTypes, ...
                'VariableNames',varNames);

    fprintf('INFO: Scanning %s\n', config.baseDir);
    fprintf(fileID, 'INFO: Scanning %s\n', config.baseDir);

    sensed_files = dir(fullfile(config.baseDir, '**', [config.file_2_name_base, '.wav']));
    fprintf('INFO: Found %d sensed files.\n', numel(sensed_files));
    fprintf(fileID, 'INFO: Found %d sensed files.\n', numel(sensed_files));

    % Case-insensitive token finder (matches path segment boundaries)
    hasTok = @(p,tok) ~isempty(regexpi(p, ['(^|[\\/_.-])' tok '([\\/_.-]|$)'], 'once'));

    for k = 1:numel(sensed_files)
        thisFolder = sensed_files(k).folder;

        % Parse Site / Day / Time from path
        site = '';
        for i = 1:numel(config.allowed_sites)
            if hasTok(thisFolder, config.allowed_sites{i})
                site = upper(config.allowed_sites{i});
                break;
            end
        end

        day = '';
        for i = 1:numel(config.allowed_days)
            if hasTok(thisFolder, config.allowed_days{i})
                day = upper(config.allowed_days{i});
                break;
            end
        end

        time = '';
        for i = 1:numel(config.allowed_times)
            if hasTok(thisFolder, config.allowed_times{i})
                time = upper(config.allowed_times{i});
                break;
            end
        end

        % Skip if any token is missing
        if any(cellfun(@isempty, {site, day, time}))
            continue;
        end

        % File paths
        f2 = fullfile(thisFolder, [config.file_2_name_base, '.wav']); % sensed
        f1 = fullfile(thisFolder, [config.file_1_name_base, '.wav']); % reference

        if ~isfile(f1)
            fprintf('WARN: Missing reference WAV for %s\n', f2);
            fprintf(fileID, 'WARN: Missing reference WAV for %s\n', f2);
            continue;
        end

        % "Trial" label = last folder name under thisFolder (robust default)
        [~, leaf] = fileparts(thisFolder);

        try
            % Compute ENF-based correlation using external pipeline
            r = proc_enf_analysis( ...
                    f1, f2, ...
                    config.nfft, config.frame_size, config.overlap_size, ...
                    config.harmonics_arr_1, config.nominal_freq_arr(2), ...
                    config.harmonics_arr_2, config.nominal_freq_arr(2), ...
                    config.trace_1_freq_est_method, config.trace_1_est_freq, config.trace_1_freq_est_spec_comb_harmonics, ...
                    config.trace_2_freq_est_method, config.trace_2_est_freq, config.trace_2_freq_est_spec_comb_harmonics, ...
                    config.trace_1_plot_title, config.trace_2_plot_title, false);

            % Add to table
            newRow = {r, string(f2), string(site), string(day), string(time)};
            tbl = [tbl; newRow];

            % Pretty log line (filenames only, no full paths)
            [~, f1name, f1ext] = fileparts(f1);
            [~, f2name, f2ext] = fileparts(f2);
            file_A_name = [f1name f1ext];
            file_B_name = [f2name f2ext];

            log_msg = sprintf(['    > Compared (%s) vs. (%s): r = %.4f  ', ...
                               '[Site=%s | Day=%s | Time=%s | Trial=%s]  ', ...
                               '{%d/%d}\n'], ...
                               file_A_name, file_B_name, r, ...
                               site, day, time, leaf, ...
                               k, numel(sensed_files));
            fprintf('%s', log_msg);
            fprintf(fileID, '%s', log_msg);

        catch ME
            fprintf('ERROR: %s\n', ME.message);
            fprintf(fileID, 'ERROR: %s\n', ME.message);
        end
    end

    fprintf('\nINFO: Batch complete. N=%d\n', height(tbl));
    fprintf(fileID, '\nINFO: Batch complete. N=%d\n', height(tbl));
end

function [r_stats, z_data, z_stats] = calculate_basic_stats(all_r)
    r_stats = calculate_descriptive_stats(all_r);
    r_bounded = all_r;
    r_bounded(r_bounded >=  1) = 0.9999999;
    r_bounded(r_bounded <= -1) = -0.9999999;
    z_data  = atanh(r_bounded);
    z_stats = calculate_descriptive_stats(z_data);
end

function stats = calculate_descriptive_stats(x)
    if ~isempty(x)
        stats.N      = numel(x);
        stats.mean   = mean(x);
        stats.median = median(x);
        stats.std    = std(x);
        stats.min    = min(x);
        stats.max    = max(x);
        stats.var    = var(x);
        stats.skew   = skewness(x);
        stats.kurt   = kurtosis(x);
        stats.sem    = stats.std / sqrt(stats.N);
        df = max(stats.N - 1, 1);
        tcrit = tinv(0.975, df);
        stats.ci95_margin = tcrit * stats.sem;
    else
        stats = struct('N',0,'mean',NaN,'median',NaN,'std',NaN,'min',NaN,'max',NaN, ...
                       'var',NaN,'skew',NaN,'kurt',NaN,'sem',NaN,'ci95_margin',NaN);
    end
end

function print_descriptive_stats_results(r_stats, z_stats, fileID)
    fprintf('\nSTEP 1A: Descriptive Results (r)\n');
    fprintf(fileID, '\nSTEP 1A: Descriptive Results (r)\n');
    print_stats_table(r_stats, fileID, 'All Correlations (r)');

    fprintf('\nSTEP 1B: Descriptive Results (z)\n');
    fprintf(fileID, '\nSTEP 1B: Descriptive Results (z)\n');
    print_stats_table(z_stats, fileID, 'All Correlations (z'')');

    fprintf('STEP 1 complete.\n');
    fprintf(fileID, 'STEP 1 complete.\n');
end

function print_stats_table(stats, fileID, label)
    fprintf('%-12s | %20s\n', 'Statistic', label);
    fprintf(fileID, '%-12s | %20s\n', 'Statistic', label);

    fprintf('%-12s | %20d\n', 'Count (N)', stats.N);
    fprintf(fileID, '%-12s | %20d\n', 'Count (N)', stats.N);

    fprintf('%-12s | %20.4f\n', 'Mean', stats.mean);
    fprintf(fileID, '%-12s | %20.4f\n', 'Mean', stats.mean);

    fprintf('%-12s | %20.4f\n', 'Median', stats.median);
    fprintf(fileID, '%-12s | %20.4f\n', 'Median', stats.median);

    fprintf('%-12s | %20.4f\n', 'Std. Dev.', stats.std);
    fprintf(fileID, '%-12s | %20.4f\n', 'Std. Dev.', stats.std);

    fprintf('%-12s | %20.4f\n', 'Variance', stats.var);
    fprintf(fileID, '%-12s | %20.4f\n', 'Variance', stats.var);

    fprintf('%-12s | %20.4f\n', 'Skewness', stats.skew);
    fprintf(fileID, '%-12s | %20.4f\n', 'Skewness', stats.skew);

    fprintf('%-12s | %20.4f\n', 'Kurtosis', stats.kurt);
    fprintf(fileID, '%-12s | %20.4f\n', 'Kurtosis', stats.kurt);

    fprintf('%-12s | %20.4f\n', 'Std. Error', stats.sem);
    fprintf(fileID, '%-12s | %20.4f\n', 'Std. Error', stats.sem);

    ci_str = sprintf('[%.4f, %.4f]', ...
        stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin);
    fprintf('%-12s | %20s\n', '95% CI', ci_str);
    fprintf(fileID, '%-12s | %20s\n', '95% CI', ci_str);

    if contains(label,"z'")
        ci_r = tanh([stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin]);
        ci_r_str = sprintf('[%.4f, %.4f]', ci_r(1), ci_r(2));
        fprintf('%-12s | %20s\n', '  (95% CI r)', ci_r_str);
        fprintf(fileID, '%-12s | %20s\n', '  (95% CI r)', ci_r_str);
    end

    fprintf('%-12s | %20.4f\n', 'Min', stats.min);
    fprintf(fileID, '%-12s | %20.4f\n', 'Min', stats.min);

    fprintf('%-12s | %20.4f\n', 'Max', stats.max);
    fprintf(fileID, '%-12s | %20.4f\n', 'Max', stats.max);
end

function plot_box_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val, 'z')
        plot_name = 'Server-Room Robustness: Fisher z-transformed Correlation Values Box Plot';
        y_label   = 'Fisher z-transform';
        suffix    = '_z';
    else
        plot_name = 'Server-Room Robustness: Pearsons Correlation Values Box Plot';
        y_label   = 'Pearson Correlation (r)';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible', 'off', 'Position', [100,100,700,600]);
    h = gca;

    if ~isempty(data)
        boxplot(h, data', 'Notch', 'on', 'Labels', {'All'});
        margin = (stats.max - stats.min) * 0.1;
        if margin == 0 || isnan(margin), margin = 0.001; end
        ylim(h, [stats.min - margin, stats.max + margin]);
        title(h, plot_name);
        ylabel(h, y_label);
        grid(h, 'on');
    else
        title(h, 'All Correlations (No Data)');
        set(h, 'XTick', [], 'YTick', []);
    end

    fn = fullfile(results_dir, ['boxplot' suffix]);
    saveas(fig, [fn '.svg']);
    fprintf('Saved boxplot%s.{svg}\n', suffix);
    fprintf(fileID, 'Saved boxplot%s.{svg}\n', suffix);
    close(fig);
end

function plot_hist_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val, 'z')
        plot_name = 'Server-Room Robustness: Fisher z-transformed Correlation Values Histogram';
        x_label   = 'Fisher z-transform';
        suffix    = '_z';
    else
        plot_name = 'Server-Room Robustness: Pearsons Correlation Values Histogram';
        x_label   = 'Pearson Correlation (r)';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible', 'off', 'Position', [100,100,700,500]);
    h = gca;

    if ~isempty(data)
        histogram(h, data, 'BinMethod', 'auto');
        margin = (stats.max - stats.min) * 0.1;
        if margin == 0 || isnan(margin), margin = 0.001; end
        xlim(h, [stats.min - margin, stats.max + margin]);
    end

    title(h, [plot_name ' (Zoomed)']);
    xlabel(h, x_label);
    ylabel(h, 'Count');
    grid(h, 'on');

    fn = fullfile(results_dir, ['hist' suffix]);
    saveas(fig, [fn '.svg']);
    fprintf('Saved hist%s.{svg}\n', suffix);
    fprintf(fileID, 'Saved hist%s.{svg}\n', suffix);
    close(fig);
end

function plot_z_score_distributions(z_data, z_stats, results_dir, fileID)
    z_mean = z_stats.mean;
    z_std  = z_stats.std;
    if z_std == 0 || isnan(z_std), z_std = 1e-6; end

    SE    = z_stats.sem;
    df    = max(z_stats.N - 1, 1);
    tcrit = tinv(0.975, df);
    CI    = [z_mean - tcrit*SE, z_mean + tcrit*SE];

    x = linspace(z_mean - 4*z_std, z_mean + 4*z_std, 1000);
    if ~isfinite(x(1)) || ~isfinite(x(end)) || numel(x) < 2 || x(1) == x(2)
        x = linspace(-3, 3, 1000);
    end

    pdf_data = normpdf(x, z_mean, z_std);

    fig = figure('Name','Z-Transformed Distributions', ...
                 'Visible','off', 'Position',[100,100,1000,600]);
    h = gca; hold(h, 'on');
    ph = plot(h, x, pdf_data, 'LineWidth', 2);
    yl = ylim(h);
    line(h, [z_mean z_mean], yl, 'LineStyle','--', 'LineWidth',1.5);
    line(h, [CI(1) CI(1)], yl, 'LineWidth',2);
    line(h, [CI(2) CI(2)], yl, 'LineWidth',2);
    ylim(h, yl);
    title(h,'Server-Room Robustness: Normal Distribution of Fisher z-transformed Correlation Values');
    xlabel(h,'z');
    ylabel(h,'PDF');
    legend(h, ph, {sprintf('PDF (mu=%.2f, std=%.2f)', z_mean, z_std)}, ...
           'Location','southeastoutside');
    grid(h,'on');
    hold(h,'off');

    fn = fullfile(results_dir,'z_score_distributions');
    saveas(fig, [fn '.svg']);
    fprintf('Saved z_score_distributions.{svg}\n');
    fprintf(fileID, 'Saved z_score_distributions.{svg}\n');
    close(fig);
end

function print_raw_correlation_data_noID(all_results_table, all_correlations, fileID)
    fprintf('All correlations calculated.\n\n');
    fprintf(fileID, 'All correlations calculated.\n\n');

    fprintf('\nFull Experiment Summary\n');
    fprintf(fileID, '\nFull Experiment Summary\n');

    try
        unique_paths = unique(all_results_table.FilePath);
        fprintf('  Total Unique Traces: %d\n', numel(unique_paths));
        fprintf(fileID,'  Total Unique Traces: %d\n', numel(unique_paths));
        total_pairs = height(all_results_table);
        fprintf('  Total Pairwise Comparisons: %d\n', total_pairs);
        fprintf(fileID,'  Total Pairwise Comparisons: %d\n', total_pairs);
    catch
    end

    fprintf('\nAll Correlations (N=%d)\n', numel(all_correlations));
    fprintf(fileID, '\nAll Correlations (N=%d)\n', numel(all_correlations));

    for i = 1:numel(all_correlations)
        fprintf('%f  ', all_correlations(i));
        fprintf(fileID, '%f  ', all_correlations(i));
        if mod(i,5) == 0
            fprintf('\n');
            fprintf(fileID, '\n');
        end
    end
    fprintf('\n');
    fprintf(fileID, '\n');

    fprintf('\nFull Results Table\n');
    fprintf(fileID, '\nFull Results Table\n');
    try
        T_str = evalc('disp(all_results_table)');
        fprintf('%s\n', T_str);
        fprintf(fileID, '%s\n', T_str);
    catch
    end
end

function debug_dump_levels(tbl, config, fileID)
    fprintf('\n[DEBUG] Levels present in tbl vs expected\n');
    fprintf(fileID, '\n[DEBUG] Levels present in tbl vs expected\n');

    sites = unique(tbl.Site);
    days  = unique(tbl.Day);
    times = unique(tbl.Time);

    fprintf('  Sites found: %s\n', strjoin(string(sites), ', '));
    fprintf('  Days  found: %s\n', strjoin(string(days),  ', '));
    fprintf('  Times found: %s\n', strjoin(string(times), ', '));

    fprintf(fileID, '  Sites found: %s\n', strjoin(string(sites), ', '));
    fprintf(fileID, '  Days  found: %s\n', strjoin(string(days),  ', '));
    fprintf(fileID, '  Times found: %s\n', strjoin(string(times), ', '));
end

function run_oneway_rm_day(tbl, config, fileID)
    T = tbl;
    T.Site = upper(string(T.Site));
    T.Day  = upper(string(T.Day));
    T.Time = upper(string(T.Time));

    keep_days  = intersect(upper(string(config.allowed_days)),  unique(T.Day));
    keep_sites = intersect(upper(string(config.allowed_sites)), unique(T.Site));
    keep_times = intersect(upper(string(config.allowed_times)), unique(T.Time));

    T = T(ismember(T.Day, keep_days) & ...
          ismember(T.Site, keep_sites) & ...
          ismember(T.Time, keep_times), :);

    if height(T) == 0
        fprintf('  [DAY] No data after filtering.\n');
        fprintf(fileID, '  [DAY] No data after filtering.\n');
        return;
    end

    rz = max(min(T.Correlation,0.9999999), -0.9999999);
    T.z = atanh(rz);

    % Mean z per Time x Day (averaging over sites)
    G = groupsummary(T, {'Time','Day'}, 'mean', 'z');
    W = unstack(G, 'mean_z', 'Day', 'GroupingVariables', 'Time'); % MON, WED

    if ~ismember('MON', W.Properties.VariableNames), W.MON = NaN(height(W),1); end
    if ~ismember('WED', W.Properties.VariableNames), W.WED = NaN(height(W),1); end

    W = W(:, {'Time','MON','WED'});
    W = rmmissing(W, 'DataVariables', {'MON','WED'});

    nSubj = height(W);
    if nSubj < 2
        fprintf('  [DAY] Not enough subjects for RM-ANOVA.\n');
        fprintf(fileID, '  [DAY] Not enough subjects for RM-ANOVA.\n');
        debug_dump_levels(T, config, fileID);
        return;
    end

    if config.SCR_CFG_LIST_SUBJECT_KEYS
        fprintf('[DAY] Subjects (Time) used:\n');
        fprintf(fileID, '[DAY] Subjects (Time) used:\n');
        disp(W.Time);
        fprintf(fileID, '%s\n', evalc('disp(W.Time)'));
    end

    Y = W{:, {'MON','WED'}};
    tY = array2table(Y, 'VariableNames', {'MON','WED'});
    within = table(categorical({'MON';'WED'}), 'VariableNames', {'Day'});

    % Assumption checks
    try
        y  = Y(:);
        id = repelem((1:nSubj)', 2);
        mu_id = grpstats(y, id, 'mean');
        resid = y - mu_id(id);
        [~, pNorm] = adtest((resid-mean(resid))/std(resid));
        g = [ones(nSubj,1); 2*ones(nSubj,1)]; % MON=1, WED=2
        pLev = vartestn(y, g, 'TestType','LeveneAbsolute','Display','off');
        fprintf('  [DAY] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
        fprintf(fileID, '  [DAY] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
    catch
    end

    % 1-way RM-ANOVA
    rm  = fitrm(tY, 'MON-WED ~ 1', 'WithinDesign', within);
    ran = ranova(rm, 'WithinModel','Day');

    rn  = string(ran.Properties.RowNames);
    idx = find(rn=="(Intercept):Day",1);
    if isempty(idx), idx = find(contains(rn,'Day'),1); end
    errRow = find(contains(rn,'Error(Day)'),1);

    df1  = ran.DF(idx);
    df2  = ran.DF(errRow);
    Fval = ran.F(idx);
    p    = ran.pValue(idx);

    if ismember('pValueGG', ran.Properties.VariableNames) && ~isnan(ran.pValueGG(idx))
        pGG = ran.pValueGG(idx);
    else
        pGG = p;
    end

    SS_eff = ran.SumSq(idx);
    SS_err = ran.SumSq(errRow);
    eta_p2 = SS_eff / (SS_eff + SS_err);

    % Print compact ANOVA summary table (ASCII)
    fprintf('    RM-ANOVA Summary (Day)\n');
    fprintf('    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf('    Day      | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    fprintf(fileID, '    RM-ANOVA Summary (Day)\n');
    fprintf(fileID, '    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf(fileID, '    Day      | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    % Level means + CI on r-scale and dot-whisker plot
    levels = {'MON','WED'};
    [mu_r, lo_r, hi_r] = print_level_means_ci_from_Y('Day', levels, Y, fileID);
    plot_dotwhisker_r(mu_r, lo_r, hi_r, levels, config.results_dir, 'day', ...
                      'Estimates with 95% CI — DAY (MON vs WED)');
end

function run_oneway_rm_site(tbl, config, fileID)
    T = tbl;
    T.Site = upper(string(T.Site));
    T.Day  = upper(string(T.Day));
    T.Time = upper(string(T.Time));

    keep_days  = intersect(upper(string(config.allowed_days)),  unique(T.Day));
    keep_sites = intersect(upper(string(config.allowed_sites)), unique(T.Site));
    keep_times = intersect(upper(string(config.allowed_times)), unique(T.Time));

    T = T(ismember(T.Day, keep_days) & ...
          ismember(T.Site, keep_sites) & ...
          ismember(T.Time, keep_times), :);

    if height(T) == 0
        fprintf('  [SITE] No data after filtering.\n');
        fprintf(fileID, '  [SITE] No data after filtering.\n');
        return;
    end

    rz = max(min(T.Correlation,0.9999999), -0.9999999);
    T.z = atanh(rz);

    % Mean z per Time x Site (averaging over days)
    G = groupsummary(T, {'Time','Site'}, 'mean', 'z');
    W = unstack(G, 'mean_z', 'Site', 'GroupingVariables', 'Time'); % SBOX, SPSU

    if ~ismember('SBOX', W.Properties.VariableNames), W.SBOX = NaN(height(W),1); end
    if ~ismember('SPSU', W.Properties.VariableNames), W.SPSU = NaN(height(W),1); end

    W = W(:, {'Time','SBOX','SPSU'});
    W = rmmissing(W, 'DataVariables', {'SBOX','SPSU'});

    nSubj = height(W);
    if nSubj < 2
        fprintf('  [SITE] Not enough subjects for RM-ANOVA.\n');
        fprintf(fileID, '  [SITE] Not enough subjects for RM-ANOVA.\n');
        debug_dump_levels(T, config, fileID);
        return;
    end

    if config.SCR_CFG_LIST_SUBJECT_KEYS
        fprintf('[SITE] Subjects (Time) used:\n');
        fprintf(fileID, '[SITE] Subjects (Time) used:\n');
        disp(W.Time);
        fprintf(fileID, '%s\n', evalc('disp(W.Time)'));
    end

    Y = W{:, {'SBOX','SPSU'}};
    tY = array2table(Y, 'VariableNames', {'SBOX','SPSU'});
    within = table(categorical({'SBOX';'SPSU'}), 'VariableNames', {'Site'});

    % Assumption checks
    try
        y  = Y(:);
        id = repelem((1:nSubj)', 2);
        mu_id = grpstats(y, id, 'mean');
        resid = y - mu_id(id);
        [~, pNorm] = adtest((resid-mean(resid))/std(resid));
        g = [ones(nSubj,1); 2*ones(nSubj,1)]; % SBOX=1, SPSU=2
        pLev = vartestn(y, g, 'TestType','LeveneAbsolute','Display','off');
        fprintf('  [SITE] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
        fprintf(fileID, '  [SITE] Assumptions: normality p=%.4f; Levene p=%.4f\n', pNorm, pLev);
    catch
    end

    % 1-way RM-ANOVA
    rm  = fitrm(tY, 'SBOX-SPSU ~ 1', 'WithinDesign', within);
    ran = ranova(rm, 'WithinModel', 'Site');

    rn  = string(ran.Properties.RowNames);
    idx = find(rn=="(Intercept):Site",1);
    if isempty(idx), idx = find(contains(rn,'Site'),1); end
    errRow = find(contains(rn,'Error(Site)'),1);

    df1  = ran.DF(idx);
    df2  = ran.DF(errRow);
    Fval = ran.F(idx);
    p    = ran.pValue(idx);

    if ismember('pValueGG', ran.Properties.VariableNames) && ~isnan(ran.pValueGG(idx))
        pGG = ran.pValueGG(idx);
    else
        pGG = p;
    end

    SS_eff = ran.SumSq(idx);
    SS_err = ran.SumSq(errRow);
    eta_p2 = SS_eff / (SS_eff + SS_err);

    % Print compact ANOVA summary table (ASCII)
    fprintf('    RM-ANOVA Summary (Site)\n');
    fprintf('    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf('    Site     | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    fprintf(fileID, '    RM-ANOVA Summary (Site)\n');
    fprintf(fileID, '    Effect   |  df1 |  df2 |        F |          p |      p(GG) |    eta_p^2\n');
    fprintf(fileID, '    Site     | %4d | %4d | %8.3f | %10.4g | %10.4g | %10.3f\n', ...
            df1, df2, Fval, p, pGG, eta_p2);

    % Level means + CI on r-scale and dot-whisker plot
    levels = {'SBOX','SPSU'};
    [mu_r, lo_r, hi_r] = print_level_means_ci_from_Y('Site', levels, Y, fileID);
    plot_dotwhisker_r(mu_r, lo_r, hi_r, levels, config.results_dir, 'site', ...
                      'Estimates with 95% CI — SITE (SBOX vs SPSU)');
end

function [mu_r, lo_r, hi_r] = print_level_means_ci_from_Y(tag, names, Y, fileID)
    % Y: nSubj x 2, z-scale values for the two levels

    n    = size(Y,1);
    mu_z = mean(Y,1);
    sd_z = std(Y,0,1);
    se_z = sd_z / max(sqrt(n),1);
    df   = max(n-1,1);
    tcrit = tinv(0.975, df);

    lo_z = mu_z - tcrit*se_z;
    hi_z = mu_z + tcrit*se_z;

    mu_r = tanh(mu_z);
    lo_r = tanh(lo_z);
    hi_r = tanh(hi_z);

    fprintf('  [%s] Level means (z''): %s=%.4f [%.4f, %.4f],  %s=%.4f [%.4f, %.4f]\n', ...
        tag, names{1}, mu_z(1), lo_z(1), hi_z(1), ...
             names{2}, mu_z(2), lo_z(2), hi_z(2));

    fprintf('  [%s] Level 95%% CI (r): %s=[%.4f,%.4f],  %s=[%.4f,%.4f]\n', ...
        tag, names{1}, lo_r(1), hi_r(1), ...
             names{2}, lo_r(2), hi_r(2));

    if fileID > 1
        fprintf(fileID, '  [%s] Level means (z''): %s=%.4f [%.4f, %.4f],  %s=%.4f [%.4f, %.4f]\n', ...
            tag, names{1}, mu_z(1), lo_z(1), hi_z(1), ...
                 names{2}, mu_z(2), lo_z(2), hi_z(2));
        fprintf(fileID, '  [%s] Level 95%% CI (r): %s=[%.4f,%.4f],  %s=[%.4f,%.4f]\n', ...
            tag, names{1}, lo_r(1), hi_r(1), ...
                 names{2}, lo_r(2), hi_r(2));
    end
end

function out = tern(c, t, f)
    if c
        out = t;
    else
        out = f;
    end
end

function plot_dotwhisker_r(mu_r, ci_lo_r, ci_hi_r, levels, results_dir, tag, title_str)
    % Dot-and-whisker plot for means ±95% CIs on r-scale.

    if nargin < 7 || isempty(title_str)
        title_str = ['Estimates with 95% CI — ' upper(tag)];
    end

    mu_r    = mu_r(:)';
    ci_lo_r = ci_lo_r(:)';
    ci_hi_r = ci_hi_r(:)';
    K       = numel(mu_r);
    x       = 1:K;

    % Simple color scheme
    mainColor    = [63  81 181] / 255;
    dotFaceColor = [26  35 126] / 255;
    dotEdgeColor = [1 1 1];

    fig = figure('Name',['DotWhisker_' tag], ...
                 'Visible','off', ...
                 'Position',[120,120,820,560]);
    ax = gca;
    hold(ax,'on');

    for i = 1:K
        line(ax, [x(i) x(i)], [ci_lo_r(i) ci_hi_r(i)], ...
            'Color', mainColor, ...
            'LineWidth', 2.0);
        plot(ax, x(i), mu_r(i), 'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', dotFaceColor, ...
            'MarkerEdgeColor', dotEdgeColor, ...
            'LineWidth', 1.2);
    end

    set(ax, 'XLim', [0.5 K+0.5], ...
            'XTick', x, ...
            'XTickLabel', levels, ...
            'FontName', 'Times', ...
            'FontSize', 10, ...
            'Box', 'off', ...
            'Layer', 'top', ...
            'YGrid', 'on', ...
            'LineWidth', 0.8);

    ylabel(ax, 'Pearson r', 'FontSize', 10);

    vals = [ci_lo_r(:); ci_hi_r(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        vals = [0.95 1.0];
    end
    ymin = max(0, min(vals) - 0.001);
    ymax = min(1.0, max(vals) + 0.001);
    if ymin >= ymax
        ymin = max(0, ymin - 0.01);
        ymax = min(1.0, ymax + 0.01);
    end
    ylim(ax, [ymin ymax]);

    step = 0.001;
    yticks = ceil(ymin/step)*step : step : floor(ymax/step)*step;
    if numel(yticks) < 2
        yticks = [ymin ymax];
    end
    set(ax, 'YTick', yticks);

    title(ax, title_str, 'FontWeight','normal', 'FontSize',10);
    grid(ax,'on');
    hold(ax,'off');

    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end

    fn = fullfile(results_dir, ['dotwhisker_' tag]);
    saveas(fig, [fn '.svg']);
    close(fig);

    fprintf('    Saved dotwhisker_%s.{svg}\n', tag);
end
