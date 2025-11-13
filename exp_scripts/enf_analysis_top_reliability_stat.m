%% ENF Temporal Reliability (Experiment 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   Quantifies how stable the ENF correlation between a sensed ambient-EM trace
%   and a mains reference is across Time-of-Day, Day-of-Week, and Week at a fixed
%   grid location. Produces summary statistics and model-based mean ±95% CI
%   (from Fisher z') using RM-ANOVA (time, day) and a paired t-test (week).
%
%ANALYSIS OVERVIEW
% 1) INPUTS
%    • Root experiment directory containing WK01, WK02, … with subfolders
%      Day ? {WED, THU, SAT, SUN}, TimeOfDay ? {EMRN, MORN, AFTN, EVEN},
%      and trial folders {T01…T05}.
%    • For each condition, two time-aligned WAV files:
%        - mains_pow_trace_ac.wav     (ground-truth mains reference)
%        - fpga_em_trace_dc.wav       (ambient EM sensed near the FPGA board)
%
% 2) ENF EXTRACTION & CORRELATION
%    • STFT spectrograms (configurable nfft/frame/overlap) around 50/60-Hz harmonics.
%    • Weighted-harmonic ENF estimator per trace.
%    • Pearson correlation r between sensed and reference ENF per file pair.
%    • Fisher transform z = atanh(r) for parametric inference.
%
% 3) STATISTICS
%    • Descriptives on r and z (mean, SD, SEM, 95% CI).
%    • Repeated-measures ANOVA on z for:
%        - TimeOfDay (EMRN, MORN, AFTN, EVEN)
%        - DayOfWeek (WED, THU, SAT, SUN)
%      Subjects = crossed (Week × other factor); assumption checks (AD normality,
%      Levene homogeneity); partial ?_p²; Greenhouse–Geisser where applicable.
%    • Week-to-Week comparison via paired t-test on z (matched Day×TimeOfDay).
%    • Holm–Bonferroni across the two RM tests.
%    • Back-transform model estimates and CIs to r for interpretability.
%
% 4) VISUALIZATION & REPORTING
%    • Dot-and-whisker plots of mean r ±95% CI for TimeOfDay, DayOfWeek, and Week.
%    • Optional histograms and z-distribution diagnostics.
%    • Optional listing of low-correlation cases below a user-set threshold.
%
% 5) OUTPUTS
%    • Logs ? exp_logs/
%    • Figures/Tables ? exp_results/enf_analysis_top_reliability_stat/
%
% Usage:
%   Set config.baseDir to the TREND root (WK01/WK02/… present) and ensure
%   proc_enf_analysis.p is on the MATLAB path. Run this script to generate the
%   full reproducible analysis and figures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% 1) Configuration
config = struct();

%Main Path Configuration
config.baseDir = '../exp_inputs/TREND/';

%File Name Configuration
config.file_1_name_base = 'mains_pow_trace_ac'; % Reference trace
config.file_2_name_base = 'fpga_em_trace_dc';   % Sensed trace

%Feature Flags
config.SHOW_RAW_CORRELATIONS      = false;   % print raw r values & full table
config.SHOW_DESCRIPTIVE_STATS     = true;    % descriptive stats + plots
config.RUN_INFERENTIAL_TESTS      = true;    % RM-ANOVA + week comparison
config.APPLY_HOLM_CORRECTION      = true;    % Holm-Bonferroni across TimeOfDay & Day
config.LIST_RM_SUBJECTS           = false;   % print subjects used in RM-ANOVA
config.SAVE_QQ_PLOTS              = false;   % save QQ plots for residuals/diffs
config.PLOT_COMBINED_INFERENTIAL  = true;    % single figure with three dot-whisker panels

%Leaf folder filter (limit to specific Txx folders)
config.FILTER_BY_LEAF_FOLDER      = true;
config.allowed_leaf_folders       = {'T01','T02','T03','T04','T05'};

%Spectrogram / ENF extraction parameters
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

%% 2) Setup: logging and colormap
[fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders();
config.results_dir = results_dir;

fprintf('Starting ENF Temporal Reliability Analysis...\n');
fprintf(fileID, 'Starting ENF Temporal Reliability Analysis. Log: %s\n\n', log_filename);

%% 3) Main analysis: compute correlations
results_tbl      = build_results_table(config, fileID);
all_correlations = results_tbl.Correlation;

fprintf('Analysis complete. N=%d correlations.\n\n', numel(all_correlations));
fprintf(fileID, 'Analysis complete. N=%d correlations.\n\n', numel(all_correlations));

%% 4) Raw results
if config.SHOW_RAW_CORRELATIONS
    print_raw_correlations(results_tbl, all_correlations, fileID);
end

%% 5) Descriptive statistics & plots
[r_stats, z_data, z_stats] = calculate_basic_stats(all_correlations);

if config.SHOW_DESCRIPTIVE_STATS
    print_descriptive_stats_results(r_stats, z_stats, fileID);
    try
        plot_box_simple(all_correlations, r_stats, results_dir, fileID, 'r');
        plot_hist_simple(all_correlations, r_stats, results_dir, fileID, 'r');
        plot_box_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_hist_simple(z_data, z_stats, results_dir, fileID, 'z');
        plot_z_score_distributions(z_data, z_stats, results_dir, fileID);
    catch MEplots
        fprintf('Plot error: %s\n', MEplots.message);
        fprintf(fileID, 'Plot error: %s\n', MEplots.message);
    end
else
    fprintf('\nSTEP 1 (Descriptive) skipped by config.\n');
    fprintf(fileID, '\nSTEP 1 (Descriptive) skipped by config.\n');
end

%% 6) Inferential statistics (TimeOfDay, Day, Week)
if config.RUN_INFERENTIAL_TESTS
    fprintf('\n\nSTEP 2: Inferential Tests (within-condition repeated measures)\n');
    fprintf(fileID, '\n\nSTEP 2: Inferential Tests (within-condition repeated measures)\n');

    try
        % Time-of-day RM-ANOVA (z-scale; report on r-scale)
        out_time = run_rm_anova(results_tbl, ...
            {'Week','Day'}, 'TimeOfDay', ...
            {'EMRN','MORN','AFTN','EVEN'}, fileID, config);

        % Day-of-week RM-ANOVA
        out_day = run_rm_anova(results_tbl, ...
            {'Week','TimeOfDay'}, 'Day', ...
            {'WED','THU','SAT','SUN'}, fileID, config);

        % Week-to-week paired comparison
        out_week = run_week_paired(results_tbl, fileID, config);

        % Holm-Bonferroni across the two RM tests (GG p-values)
        if config.APPLY_HOLM_CORRECTION
            apply_holm_two_tests(out_time, out_day, fileID);
        end

        % Combined figure with all three dot-whisker panels
        if config.PLOT_COMBINED_INFERENTIAL
            try
                plot_combined_inferential_boxplots(out_time, out_day, out_week, results_dir);
            catch MEc
                fprintf('  [WARN] Combined inferential plot failed: %s\n', MEc.message);
                fprintf(fileID, '  [WARN] Combined inferential plot failed: %s\n', MEc.message);
            end
        end

    catch MEopt
        fprintf('  [ERROR] Inferential tests failed: %s\n', MEopt.message);
        fprintf(fileID, '  [ERROR] Inferential tests failed: %s\n', MEopt.message);
    end
else
    fprintf('\n\nSTEP 2 skipped by config.\n');
    fprintf(fileID, '\n\nSTEP 2 skipped by config.\n');
end

fprintf('\n\nAll analysis complete. Log saved.\n');
fprintf(fileID, '\n\nAll analysis complete.\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedures for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fileID, results_dir, log_filename, cleanupObj] = setup_logging_and_folders()
    log_dir     = 'exp_logs';
    results_dir = 'exp_results/enf_analysis_top_reliability_stat';
    
    if ~exist(log_dir, 'dir'), mkdir(log_dir); end
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end

    [~, script_name] = fileparts(mfilename('fullpath'));
    if isempty(script_name)
        script_name = 'enf_analysis_top_reliability_stat';
    end

    log_filename = fullfile(log_dir, [script_name, '_log.txt']);
    if exist(log_filename, 'file'), delete(log_filename); end

    fileID = fopen(log_filename, 'w');
    if fileID == -1, error('Could not open %s', log_filename); end

    cleanupObj = onCleanup(@() fclose(fileID));
end

function tbl = build_results_table(config, fileID)
    % Scan for sensed files, match with mains, run ENF correlation.
    varTypes = {'double','string','string','string','string','string'};
    varNames = {'Correlation','FilePath','DayType','TimeOfDay','Week','Day'};
    tbl = table('Size',[0, numel(varNames)], ...
                'VariableTypes',varTypes, 'VariableNames',varNames);

    fprintf('INFO: Searching for trace pairs in: %s\n', config.baseDir);
    fprintf(fileID, 'INFO: Searching for trace pairs in: %s\n', config.baseDir);

    sensed_files = dir(fullfile(config.baseDir, '**', [config.file_2_name_base, '.wav']));
    fprintf('INFO: Found %d sensed files.\n', numel(sensed_files));
    fprintf(fileID, 'INFO: Found %d sensed files.\n', numel(sensed_files));

    hasToken = @(p, tok) ~isempty(regexp(p, ['(^|[\\/_.-])' tok '([\\/_.-]|$)'], 'once'));

    for k = 1:numel(sensed_files)
        currentFolder = sensed_files(k).folder;

        % Optional leaf folder filter
        if config.FILTER_BY_LEAF_FOLDER
            [~, leaf] = fileparts(currentFolder);
            if ~ismember(leaf, config.allowed_leaf_folders)
                fprintf('\n--- Skip %d/%d (Leaf: %s not in allowed list)\n', ...
                        k, numel(sensed_files), leaf);
                fprintf(fileID, '\n--- Skip %d/%d (Leaf: %s not in allowed list)\n', ...
                        k, numel(sensed_files), leaf);
                continue;
            end
        end

        f2 = fullfile(currentFolder, [config.file_2_name_base, '.wav']); % sensed
        f1 = fullfile(currentFolder, [config.file_1_name_base, '.wav']); % reference

        if ~isfile(f1)
            fprintf('WARN: Missing reference WAV at %s\n', currentFolder);
            fprintf(fileID, 'WARN: Missing reference WAV at %s\n', currentFolder);
            continue;
        end

        % Parse labels from folder path
        dayType = 'NA';
        if     hasToken(currentFolder, 'WDAY'), dayType = 'WDAY';
        elseif hasToken(currentFolder, 'WEND'), dayType = 'WEND'; end

        day = 'NA';
        if     hasToken(currentFolder, 'WED'), day = 'WED';
        elseif hasToken(currentFolder,'THU') || hasToken(currentFolder,'THUR'), day = 'THU';
        elseif hasToken(currentFolder, 'SAT'), day = 'SAT';
        elseif hasToken(currentFolder, 'SUN'), day = 'SUN'; end

        timeOfDay = 'NA';
        if     hasToken(currentFolder, 'EMRN'), timeOfDay = 'EMRN';
        elseif hasToken(currentFolder, 'MORN'),  timeOfDay = 'MORN';
        elseif hasToken(currentFolder, 'AFTN'),  timeOfDay = 'AFTN';
        elseif hasToken(currentFolder, 'EVEN'),  timeOfDay = 'EVEN'; end

        week = 'NA';
        if     hasToken(currentFolder, 'WK01'), week = 'WK01';
        elseif hasToken(currentFolder, 'WK02'), week = 'WK02'; end

        try
            r = proc_enf_analysis( ...
                    f1, f2, ...
                    config.nfft, config.frame_size, config.overlap_size, ...
                    config.harmonics_arr_1, config.nominal_freq_arr(2), ...
                    config.harmonics_arr_2, config.nominal_freq_arr(2), ...
                    config.trace_1_freq_est_method, config.trace_1_est_freq, config.trace_1_freq_est_spec_comb_harmonics, ...
                    config.trace_2_freq_est_method, config.trace_2_est_freq, config.trace_2_freq_est_spec_comb_harmonics, ...
                    config.trace_1_plot_title, config.trace_2_plot_title, false);

            newRow = {r, string(f2), string(dayType), string(timeOfDay), string(week), string(day)};
            tbl    = [tbl; newRow];
            
            % Build nice log line with just the file basenames
            [~, f1name, f1ext] = fileparts(f1);
            [~, f2name, f2ext] = fileparts(f2);
            file_A_name = [f1name f1ext];
            file_B_name = [f2name f2ext];

            log_msg = sprintf(['    > Compared (%s) vs. (%s): r = %.4f  ', ...
                               '[Week=%s | Day=%s | Time=%s | Trial=%s]  ', ...
                               '{%d/%d}\n'], ...
                               file_A_name, file_B_name, r, ...
                               week, day, timeOfDay, leaf, ...
                               k, numel(sensed_files));
            fprintf('%s', log_msg);
            fprintf(fileID, '%s', log_msg);

            fprintf('INFO: %d/%d corr = %.6f\n', k, numel(sensed_files), r);
            fprintf(fileID, 'INFO: %d/%d corr = %.6f\n', k, numel(sensed_files), r);
        catch ME
            fprintf('ERROR: Analysis failed for pair %d: %s\n', k, ME.message);
            fprintf(fileID, 'ERROR: Analysis failed for pair %d: %s\n', k, ME.message);
        end
    end

    fprintf('\nINFO: Batch complete. N=%d\n', height(tbl));
    fprintf(fileID, '\nINFO: Batch complete. N=%d\n', height(tbl));
end

function [r_stats, z_data, z_stats] = calculate_basic_stats(all_correlations)
    r_stats = calculate_descriptive_stats(all_correlations);

    r_bounded = all_correlations;
    r_bounded(r_bounded >= 1.0)  = 0.9999999;
    r_bounded(r_bounded <= -1.0) = -0.9999999;

    z_data  = atanh(r_bounded);
    z_stats = calculate_descriptive_stats(z_data);
end

function stats = calculate_descriptive_stats(x)
    if isempty(x)
        stats = struct('N',0,'mean',NaN,'median',NaN,'std',NaN,'min',NaN,'max',NaN, ...
                       'var',NaN,'skew',NaN,'kurt',NaN,'sem',NaN,'ci95_margin',NaN);
        return;
    end
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
    df           = max(stats.N - 1, 1);
    tcrit        = tinv(0.975, df);
    stats.ci95_margin = tcrit * stats.sem;
end

function print_descriptive_stats_results(r_stats, z_stats, fileID)
    % STEP 1A: r-scale
    fprintf('\nSTEP 1A: Descriptive Results (r)\n');
    fprintf(fileID, '\nSTEP 1A: Descriptive Results (r)\n');
    print_stats_table(r_stats, fileID, 'All Correlations (r)');

    % STEP 1B: Fisher z
    fprintf('\nSTEP 1B: Descriptive Results (z)\n');
    fprintf(fileID, '\nSTEP 1B: Descriptive Results (z)\n');

    % back-transformed mean and CI on r-scale
    z_mean        = z_stats.mean;
    mean_r_from_z = tanh(z_mean);
    ci_z          = [z_stats.mean - z_stats.ci95_margin, z_stats.mean + z_stats.ci95_margin];
    ci_r_from_z   = tanh(ci_z);

    print_stats_table(z_stats, fileID, 'All Correlations (z'')');
    fprintf('Mean (r from z-mean)   | %20.4f\n', mean_r_from_z);
    fprintf('95%% CI (r from z)      | [%0.4f, %0.4f]\n', ci_r_from_z(1), ci_r_from_z(2));
    fprintf(fileID, 'Mean (r from z-mean)   | %20.4f\n', mean_r_from_z);
    fprintf(fileID, '95%% CI (r from z)      | [%0.4f, %0.4f]\n', ci_r_from_z(1), ci_r_from_z(2));

    fprintf('STEP 1 complete.\n');
    fprintf(fileID, 'STEP 1 complete.\n');
end

function print_stats_table(stats, fileID, label)
    fprintf('%-20s | %20s\n', 'Statistic', label);
    fprintf(fileID, '%-20s | %20s\n', 'Statistic', label);

    fprintf('%-20s | %20d\n', 'Count (N)', stats.N);
    fprintf('%-20s | %20.4f\n', 'Mean', stats.mean);
    fprintf('%-20s | %20.4f\n', 'Median', stats.median);
    fprintf('%-20s | %20.4f\n', 'Std. Dev.', stats.std);
    fprintf('%-20s | %20.4f\n', 'Variance', stats.var);
    fprintf('%-20s | %20.4f\n', 'Skewness', stats.skew);
    fprintf('%-20s | %20.4f\n', 'Kurtosis', stats.kurt);
    fprintf('%-20s | %20.4f\n', 'Std. Error', stats.sem);
    fprintf('%-20s | [%0.4f, %0.4f]\n', '95% CI', ...
            stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin);
    fprintf('%-20s | %20.4f\n', 'Min', stats.min);
    fprintf('%-20s | %20.4f\n', 'Max', stats.max);

    fprintf(fileID, '%-20s | %20d\n', 'Count (N)', stats.N);
    fprintf(fileID, '%-20s | %20.4f\n', 'Mean', stats.mean);
    fprintf(fileID, '%-20s | %20.4f\n', 'Median', stats.median);
    fprintf(fileID, '%-20s | %20.4f\n', 'Std. Dev.', stats.std);
    fprintf(fileID, '%-20s | %20.4f\n', 'Variance', stats.var);
    fprintf(fileID, '%-20s | %20.4f\n', 'Skewness', stats.skew);
    fprintf(fileID, '%-20s | %20.4f\n', 'Kurtosis', stats.kurt);
    fprintf(fileID, '%-20s | %20.4f\n', 'Std. Error', stats.sem);
    fprintf(fileID, '%-20s | [%0.4f, %0.4f]\n', '95% CI', ...
            stats.mean - stats.ci95_margin, stats.mean + stats.ci95_margin);
    fprintf(fileID, '%-20s | %20.4f\n', 'Min', stats.min);
    fprintf(fileID, '%-20s | %20.4f\n', 'Max', stats.max);
end

function out = run_rm_anova(tbl, subjectGroupingVars, factorName, levels, fileID, config)
    % RM-ANOVA on Fisher z; plots factor-level dot-whisker (mean r ±95% CI).

    % Prepare Fisher z
    r_bounded = tbl.Correlation;
    r_bounded(r_bounded >= 1.0)  = 0.9999999;
    r_bounded(r_bounded <= -1.0) = -0.9999999;
    tbl.z = atanh(r_bounded);

    tbl.Day       = categorical(tbl.Day);
    tbl.Week      = categorical(tbl.Week);
    tbl.TimeOfDay = categorical(tbl.TimeOfDay);

    fprintf('\n[%s] Repeated-Measures ANOVA (subjects = %s)\n', ...
            factorName, strjoin(subjectGroupingVars,' x '));
    fprintf(fileID, '\n[%s] Repeated-Measures ANOVA (subjects = %s)\n', ...
            factorName, strjoin(subjectGroupingVars,' x '));

    % Build mean z per subject x factor level
    S = groupsummary(tbl, [subjectGroupingVars, {factorName}], 'mean', 'z');
    W = unstack(S, 'mean_z', factorName, 'GroupingVariables', subjectGroupingVars);

    % Ensure all level columns exist
    for i = 1:numel(levels)
        if ~ismember(levels{i}, W.Properties.VariableNames)
            W.(levels{i}) = nan(height(W),1);
        end
    end

    % Retain complete subjects
    W = rmmissing(W, 'DataVariables', levels);
    nSubj = height(W);

    if config.LIST_RM_SUBJECTS
        fprintf('  Subjects retained for %s = %d\n', factorName, nSubj);
        fprintf(fileID, '  Subjects retained for %s = %d\n', factorName, nSubj);
        try
            keyTbl = W(:, subjectGroupingVars);
            keyStr = evalc('disp(keyTbl)');
            fprintf('%s\n', keyStr);
            fprintf(fileID, '%s\n', keyStr);
        catch
        end
    end

    if nSubj < 2
        fprintf('  [SKIP] Not enough subjects.\n');
        fprintf(fileID, '  [SKIP] Not enough subjects.\n');
        out = struct('factor',factorName,'p_raw',NaN,'p_GG',NaN,'eta_p2',NaN, ...
                     'n',nSubj,'levels',{levels}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
        return;
    end

    Y  = W{:, levels};        % rows: subjects, cols: levels (Fisher z)
    tY = array2table(Y, 'VariableNames', levels);
    withinDesign = table(categorical(levels'), 'VariableNames', {factorName});

    % Assumption checks
    try
        y_stack = Y(:);
        g       = repelem(categorical(levels'), nSubj);
        id      = repelem((1:nSubj)', numel(levels));
        mu_id   = grpstats(y_stack, id, 'mean');
        resid   = y_stack - mu_id(id);

        xs = (resid - mean(resid)) / std(resid);
        [~, pNorm] = adtest(xs);
        pLevene = vartestn(y_stack, g, 'TestType','LeveneAbsolute','Display','off');

        fprintf('    Assumptions: normality p=%.4f (AD); Levene p=%.4f\n', pNorm, pLevene);
        fprintf(fileID, '    Assumptions: normality p=%.4f (AD); Levene p=%.4f\n', pNorm, pLevene);

        if config.SAVE_QQ_PLOTS
            fig = figure('Name', ['QQ Residuals: ' factorName], ...
                         'Visible','off', 'Position',[120,120,700,550]);
            qqplot(resid);
            title(['QQ Plot of Residuals - ' factorName]);
            xlabel('Theoretical Quantiles'); ylabel('Sample Quantiles'); grid on;
            qq_fn = fullfile(config.results_dir, ['qq_residuals_' lower(factorName)]);
            saveas(fig, [qq_fn '.png']);
            saveas(fig, [qq_fn '.svg']);
            savefig(fig, [qq_fn '.fig']);
            close(fig);
        end
    catch MEa
        fprintf('    [WARN] Assumption checks skipped: %s\n', MEa.message);
        fprintf(fileID, '    [WARN] Assumption checks skipped: %s\n', MEa.message);
    end

    try
        % Fit RM model
        formula = sprintf('%s-%s ~ 1', levels{1}, levels{end});
        rm  = fitrm(tY, formula, 'WithinDesign', withinDesign);
        ran = ranova(rm, 'WithinModel', factorName);

        % Locate factor and error rows
        rn       = string(ran.Properties.RowNames);
        termName = "(Intercept):" + string(factorName);
        idx      = find(rn == termName, 1);
        if isempty(idx)
            idx = find(contains(rn, string(factorName)), 1, 'first');
        end
        errIdx = find(contains(rn, "Error(" + string(factorName) + ")"), 1, 'first');
        if isempty(errIdx), errIdx = idx + 1; end
        if isempty(idx)
            error('Could not locate factor row for %s.', factorName);
        end

        % F, p, effect size
        df1   = ran.DF(idx);
        df2   = ran.DF(errIdx);
        Fval  = ran.F(idx);
        p_raw = ran.pValue(idx);

        p_GG = p_raw;
        if ismember('pValueGG', ran.Properties.VariableNames) && ~isnan(ran.pValueGG(idx))
            p_GG = ran.pValueGG(idx);
        end

        SS_eff = ran.SumSq(idx);
        SS_err = ran.SumSq(errIdx);
        eta_p2 = SS_eff / (SS_eff + SS_err);

        % Optional epsilon (log only)
        try
            epsTbl = epsilon(rm);
            eRow = epsTbl(strcmp(string(epsTbl.WithinEffect), factorName), :);
        catch
            eRow = [];
        end

        % RM-ANOVA summary table
        fprintf('\n    RM-ANOVA Summary (%s)\n', factorName);
        fprintf('    %-10s | %4s | %4s | %8s | %10s | %10s | %10s\n', ...
            'Effect','df1','df2','F','p','p(GG)','eta_p^2');
        fprintf('    %-10s | %4d | %4d | %8.3f | %10.4f | %10.4f | %10.3f\n', ...
            factorName, df1, df2, Fval, p_raw, p_GG, eta_p2);

        fprintf(fileID, '\n    RM-ANOVA Summary (%s)\n', factorName);
        fprintf(fileID, ...
            '    %-10s | %4s | %4s | %8s | %10s | %10s | %10s\n', ...
            'Effect','df1','df2','F','p','p(GG)','eta_p^2');
        fprintf(fileID, ...
            '    %-10s | %4d | %4d | %8.3f | %10.4f | %10.4f | %10.3f\n', ...
            factorName, df1, df2, Fval, p_raw, p_GG, eta_p2);

        if ~isempty(eRow)
            fprintf('    Epsilon (GG/HF/LB) = [%.3f / %.3f / %.3f]\n', ...
                eRow.GreenhouseGeisser, eRow.HuynhFeldt, eRow.LowerBound);
            fprintf(fileID, '    Epsilon (GG/HF/LB) = [%.3f / %.3f / %.3f]\n', ...
                eRow.GreenhouseGeisser, eRow.HuynhFeldt, eRow.LowerBound);
        end

        %-- Level means and 95% CI on z-scale--
        mu_z    = mean(Y,1);
        sd_z    = std(Y,0,1);
        se_z    = sd_z ./ sqrt(nSubj);
        dfL     = max(nSubj - 1,1);
        tcrit   = tinv(0.975, dfL);
        ci_lo_z = mu_z - tcrit .* se_z;
        ci_hi_z = mu_z + tcrit .* se_z;

        % Back-transform to r-scale
        mu_r    = tanh(mu_z);
        ci_lo_r = tanh(ci_lo_z);
        ci_hi_r = tanh(ci_hi_z);

        fprintf('    Level means (z''''): ');
        for j = 1:numel(levels), fprintf('%s=%.3f ', levels{j}, mu_z(j)); end
        fprintf('\n');

        fprintf('    Level means (r from z): ');
        for j = 1:numel(levels), fprintf('%s=%.4f ', levels{j}, mu_r(j)); end
        fprintf('\n');

        fprintf('    Level 95%% CI (r from z): ');
        for j = 1:numel(levels)
            fprintf('%s=[%.4f,%.4f] ', levels{j}, ci_lo_r(j), ci_hi_r(j));
        end
        fprintf('\n');

        fprintf(fileID, '    Level means (z''''): ');
        for j = 1:numel(levels), fprintf(fileID, '%s=%.3f ', levels{j}, mu_z(j)); end
        fprintf(fileID, '\n');

        fprintf(fileID, '    Level means (r from z): ');
        for j = 1:numel(levels), fprintf(fileID, '%s=%.4f ', levels{j}, mu_r(j)); end
        fprintf(fileID, '\n');

        fprintf(fileID, '    Level 95%% CI (r from z): ');
        for j = 1:numel(levels)
            fprintf(fileID, '%s=[%.4f,%.4f] ', levels{j}, ci_lo_r(j), ci_hi_r(j));
        end
        fprintf(fileID, '\n');

        out = struct( ...
            'factor',   factorName, ...
            'p_raw',    p_raw, ...
            'p_GG',     p_GG, ...
            'eta_p2',   eta_p2, ...
            'n',        nSubj, ...
            'levels',   {levels}, ...
            'mu_r',     mu_r, ...
            'ci_lo_r',  ci_lo_r, ...
            'ci_hi_r',  ci_hi_r);

    catch MEf
        fprintf('    [ERROR] RM-ANOVA failed for %s: %s\n', factorName, MEf.message);
        fprintf(fileID, '    [ERROR] RM-ANOVA failed for %s: %s\n', factorName, MEf.message);
        out = struct('factor',factorName,'p_raw',NaN,'p_GG',NaN,'eta_p2',NaN, ...
                     'n',nSubj,'levels',{levels}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
    end
end

function out = run_week_paired(tbl, fileID, config)
    % Paired WK01 vs WK02 comparison across matched (Day, TimeOfDay).

    r_bounded = tbl.Correlation;
    r_bounded(r_bounded >= 1.0)  = 0.9999999;
    r_bounded(r_bounded <= -1.0) = -0.9999999;
    tbl.z = atanh(r_bounded);

    tbl.Day       = categorical(tbl.Day);
    tbl.Week      = categorical(tbl.Week);
    tbl.TimeOfDay = categorical(tbl.TimeOfDay);

    fprintf('\n[Week] Paired comparison WK01 vs WK02 (matched Day x TimeOfDay)\n');
    fprintf(fileID, '\n[Week] Paired comparison WK01 vs WK02 (matched Day x TimeOfDay)\n');

    S = groupsummary(tbl, {'Day','TimeOfDay','Week'}, 'mean', 'z');
    W = unstack(S, 'mean_z', 'Week', 'GroupingVariables', {'Day','TimeOfDay'});

    if ~ismember('WK01', W.Properties.VariableNames), W.WK01 = nan(height(W),1); end
    if ~ismember('WK02', W.Properties.VariableNames), W.WK02 = nan(height(W),1); end

    W = rmmissing(W, 'DataVariables', {'WK01','WK02'});
    nPairs = height(W);

    if nPairs < 2
        fprintf('  [SKIP] Not enough matched pairs.\n');
        fprintf(fileID, '  [SKIP] Not enough matched pairs.\n');
        out = struct('factor','Week','p_t',NaN,'p_wil',NaN,'df',NaN, ...
                     'n',nPairs,'levels',{{'WK01','WK02'}}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
        return;
    end

    x = W.WK01;        % mean z per condition for WK01
    y = W.WK02;        % mean z per condition for WK02
    d = x - y;

    % Normality (AD) on differences
    try
        xs = (d - mean(d)) / std(d);
        [~, pNorm] = adtest(xs);
        fprintf('    Normality (diff) p=%.4f (AD)\n', pNorm);
        fprintf(fileID, '    Normality (diff) p=%.4f (AD)\n', pNorm);
    catch MEa
        fprintf('    [WARN] Normality check skipped: %s\n', MEa.message);
        fprintf(fileID, '    [WARN] Normality check skipped: %s\n', MEa.message);
    end

    try
        [~, p_t, ci_z, stats] = ttest(x, y);
        mean_diff_z = mean(d);
        p_wil = signrank(x, y);

        % Console
        fprintf('\n    Week Comparison Summary (z-scale)\n');
        fprintf('    %-12s | %4s | %7s | %10s | %25s\n', ...
                'Test','df','t','p','Mean diff z [95% CI]');
        fprintf('    %-12s | %4d | %7.3f | %10.4f | %8.4f [%8.4f, %8.4f]\n', ...
                'paired t', stats.df, stats.tstat, p_t, mean_diff_z, ci_z(1), ci_z(2));
        fprintf('    %-12s | %4s | %7s | %10.4f | %25s\n', ...
                'Wilcoxon','--','--', p_wil, 'n/a');

        % Log file
        fprintf(fileID, '\n    Week Comparison Summary (z-scale)\n');
        fprintf(fileID, '    %-12s | %4s | %7s | %10s | %25s\n', ...
                        'Test','df','t','p','Mean diff z [95% CI]');
        fprintf(fileID, '    %-12s | %4d | %7.3f | %10.4f | %8.4f [%8.4f, %8.4f]\n', ...
                        'paired t', stats.df, stats.tstat, p_t, mean_diff_z, ci_z(1), ci_z(2));
        fprintf(fileID, '    %-12s | %4s | %7s | %10.4f | %25s\n', ...
                        'Wilcoxon','--','--', p_wil, 'n/a');


        % Back-transformed WK means & 95% CI (for plotting)
        levels = {'WK01','WK02'};
        mu_z   = [mean(x), mean(y)];
        sd_z   = [std(x),  std(y)];
        se_z   = sd_z ./ sqrt(nPairs);
        dfL    = max(nPairs - 1, 1);
        tcrit  = tinv(0.975, dfL);
        ci_lo_z = mu_z - tcrit .* se_z;
        ci_hi_z = mu_z + tcrit .* se_z;

        mu_r    = tanh(mu_z);
        ci_lo_r = tanh(ci_lo_z);
        ci_hi_r = tanh(ci_hi_z);

        fprintf('\n    Week means (r from z): WK01=%.4f, WK02=%.4f\n', ...
            mu_r(1), mu_r(2));
        fprintf(fileID, '\n    Week means (r from z): WK01=%.4f, WK02=%.4f\n', ...
            mu_r(1), mu_r(2));

        out = struct( ...
            'factor',   'Week', ...
            'p_t',      p_t, ...
            'p_wil',    p_wil, ...
            'df',       stats.df, ...
            'n',        nPairs, ...
            'levels',   {{'WK01','WK02'}}, ...
            'mu_r',     mu_r, ...
            'ci_lo_r',  ci_lo_r, ...
            'ci_hi_r',  ci_hi_r);

    catch MEp
        fprintf('    [ERROR] Paired tests failed: %s\n', MEp.message);
        fprintf(fileID, '    [ERROR] Paired tests failed: %s\n', MEp.message);
        out = struct('factor','Week','p_t',NaN,'p_wil',NaN,'df',NaN, ...
                     'n',nPairs,'levels',{{'WK01','WK02'}}, ...
                     'mu_r',[],'ci_lo_r',[],'ci_hi_r',[]);
    end
end

function plot_combined_inferential_boxplots(out_time, out_day, out_week, results_dir)
% Combined figure with three dot-and-whisker panels:
% (1) Time-of-day, (2) Day-of-week, (3) Week-to-week.
% Each panel shows mean r and 95% CI on the r-scale.
% Uses same violet color scheme as individual plots.

    fig = figure('Name','Inferential_DotWhisker_Combined', ...
                 'Visible','off', ...
                 'Position',[100,100,1350,420]);

    violet_main = [63  81 181] / 255;
    violet_dark = [26  35 126] / 255;

    %Subplot 1: Time-of-day
    ax1 = subplot(1,3,1);
    if isfield(out_time,'mu_r') && ~isempty(out_time.mu_r)
        draw_dotwhisker_on_ax(ax1, out_time.mu_r, out_time.ci_lo_r, out_time.ci_hi_r, ...
            out_time.levels, 'Time-of-day', true, violet_main, violet_dark);
    else
        axis(ax1,'off');
        text(0.5,0.5,'Time-of-day N/A','HorizontalAlignment','center');
    end

    %Subplot 2: Day-of-week
    ax2 = subplot(1,3,2);
    if isfield(out_day,'mu_r') && ~isempty(out_day.mu_r)
        draw_dotwhisker_on_ax(ax2, out_day.mu_r, out_day.ci_lo_r, out_day.ci_hi_r, ...
            out_day.levels, 'Day-of-week', false, violet_main, violet_dark);
    else
        axis(ax2,'off');
        text(0.5,0.5,'Day-of-week N/A','HorizontalAlignment','center');
    end

    %Subplot 3: Week-to-week
    ax3 = subplot(1,3,3);
    if isfield(out_week,'mu_r') && ~isempty(out_week.mu_r)
        draw_dotwhisker_on_ax(ax3, out_week.mu_r, out_week.ci_lo_r, out_week.ci_hi_r, ...
            out_week.levels, 'Week-to-week', false, violet_main, violet_dark);
    else
        axis(ax3,'off');
        text(0.5,0.5,'Week N/A','HorizontalAlignment','center');
    end

    % Save combined figure
    fn = fullfile(results_dir, 'inferential_dotwhisker_combined');
    saveas(fig, [fn '.svg']);
    close(fig);

    fprintf('  Saved inferential_dotwhisker_combined.{png,svg,fig}\n');
end

function draw_dotwhisker_on_ax(ax, mu_r, ci_lo_r, ci_hi_r, levels, title_str, showYLabel, ...
                               violet_main, violet_dark)
    axes(ax); 
    hold(ax,'on');

    mu_r    = mu_r(:)';
    ci_lo_r = ci_lo_r(:)';
    ci_hi_r = ci_hi_r(:)';
    K       = numel(mu_r);
    x       = 1:K;

    whiskerColor = violet_main;
    dotFaceColor = violet_dark;
    dotEdgeColor = [1 1 1];

    for i = 1:K
        % vertical CI line
        line(ax, [x(i) x(i)], [ci_lo_r(i) ci_hi_r(i)], ...
            'Color', whiskerColor, ...
            'LineWidth', 2.0);
        % mean dot
        plot(ax, x(i), mu_r(i), 'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', dotFaceColor, ...
            'MarkerEdgeColor', dotEdgeColor, ...
            'LineWidth', 1.2);
    end

    if K > 1
        plot(ax, x, mu_r, '-', ...
            'Color', whiskerColor, ...
            'LineWidth', 1.0);
    end

    set(ax, 'XLim', [0.5 K+0.5], ...
            'XTick', x, ...
            'XTickLabel', levels, ...
            'FontName', 'Times', ...
            'FontSize', 9, ...
            'Box', 'off', ...
            'Layer', 'top', ...
            'YGrid', 'on', ...
            'GridColor', [0.90 0.91 0.97], ...
            'GridAlpha', 0.9, ...
            'LineWidth', 0.8);

    if showYLabel
        ylabel(ax, 'Pearson r', 'FontSize', 9);
    else
        ylabel(ax, '');
    end

    % y-limits based on CI, clipped to [0,1]
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

    title(ax, title_str, 'FontWeight','normal', 'FontSize',10);
    hold(ax,'off');
end

function apply_holm_two_tests(outA, outB, fileID)
    fprintf('\n[Holm-Bonferroni across RM tests (GG p-values)]\n');
    fprintf(fileID, '\n[Holm-Bonferroni across RM tests (GG p-values)]\n');

    pvals = [outA.p_GG, outB.p_GG];
    names = {outA.factor, outB.factor};

    [p_sorted, order] = sort(pvals);
    m = numel(pvals);
    any_sig = false;

    for i = 1:m
        idx     = order(i);
        alpha_i = 0.05 / (m - i + 1);
        decision = (p_sorted(i) < alpha_i);
        any_sig  = any_sig | decision;

        fprintf('  %s: p=%.4f vs alpha_Holm=%.4f -> %s\n', ...
            names{idx}, p_sorted(i), alpha_i, ternary(decision,'SIGNIFICANT','ns'));
        fprintf(fileID, '  %s: p=%.4f vs alpha_Holm=%.4f -> %s\n', ...
            names{idx}, p_sorted(i), alpha_i, ternary(decision,'SIGNIFICANT','ns'));

        if ~decision
            fprintf('  (Holm stops here; remaining are non-significant.)\n');
            fprintf(fileID, '  (Holm stops here; remaining are non-significant.)\n');
            break;
        end
    end

    if ~any_sig
        fprintf('  Result: No RM effect survives Holm correction.\n');
        fprintf(fileID, '  Result: No RM effect survives Holm correction.\n');
    end
end

function print_raw_correlations(results_tbl, all_correlations, fileID)
    fprintf('All correlations calculated.\n\n');
    fprintf(fileID, 'All correlations calculated.\n\n');

    fprintf('\nFull Experiment Summary\n');
    fprintf(fileID, '\nFull Experiment Summary\n');

    try
        unique_paths = unique(results_tbl.FilePath);
        fprintf('  Total Unique Traces: %d\n', numel(unique_paths));
        fprintf('  Total Pairwise Comparisons: %d\n', height(results_tbl));
        fprintf(fileID, '  Total Unique Traces: %d\n', numel(unique_paths));
        fprintf(fileID, '  Total Pairwise Comparisons: %d\n', height(results_tbl));
    catch ME
        fprintf('  Summary error: %s\n', ME.message);
        fprintf(fileID, '  Summary error: %s\n', ME.message);
    end

    fprintf('\nAll Correlations (N=%d)\n', numel(all_correlations));
    fprintf(fileID, '\nAll Correlations (N=%d)\n', numel(all_correlations));
    for i = 1:numel(all_correlations)
        fprintf('%0.6f  ', all_correlations(i));
        fprintf(fileID, '%0.6f  ', all_correlations(i));
        if mod(i,5)==0
            fprintf('\n'); fprintf(fileID,'\n');
        end
    end
    fprintf('\n'); fprintf(fileID,'\n');

    fprintf('\nFull Results Table\n');
    fprintf(fileID, '\nFull Results Table\n');
    try
        T_str = evalc('disp(results_tbl)');
        fprintf('%s\n', T_str);
        fprintf(fileID, '%s\n', T_str);
    catch ME2
        fprintf('  Table print error: %s\n', ME2.message);
        fprintf(fileID, '  Table print error: %s\n', ME2.message);
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function plot_box_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val,'z')
        plot_name = 'Temporal Reliability: Fisher z-transformed Correlation Values Box Plot';
        y_label   = 'Fisher z';
        suffix    = '_z';
    else
        plot_name = 'Temporal Reliability: Pearsons Correlation Values Box Plot';
        y_label   = 'Pearson r';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible','off', ...
                 'Position',[100,100,700,600]);
    ax = gca;

    if ~isempty(data)
        boxplot(ax, data', 'Notch','on','Labels',{'All'});
        margin = (stats.max - stats.min)*0.1;
        if margin <= 0 || isnan(margin), margin = 0.001; end
        ylim(ax, [stats.min - margin, stats.max + margin]);
        ylabel(ax, y_label);
        title(ax, plot_name);
        grid(ax,'on');
    else
        title(ax,'No Data');
        set(ax,'XTick',[],'YTick',[]);
    end

    saveas(fig, fullfile(results_dir, ['boxplot' suffix '.svg']));
    fprintf('Saved boxplot%s.{svg}\n', suffix);
    fprintf(fileID,'Saved boxplot%s.{svg}\n', suffix);
    close(fig);
end

function plot_hist_simple(data, stats, results_dir, fileID, which_val)
    if strcmp(which_val,'z')
        plot_name = 'Temporal Reliability: Fisher z-transformed Correlation Values Histogram';
        x_label   = 'Fisher z';
        suffix    = '_z';
    else
        plot_name = 'Temporal Reliability: Pearsons Correlation Values Histogram';
        x_label   = 'Pearson r';
        suffix    = '_r';
    end

    fig = figure('Name', plot_name, 'Visible','off', ...
                 'Position',[100,100,700,500]);
    ax = gca;

    if ~isempty(data)
        histogram(ax, data, 'BinMethod','auto');
        margin = (stats.max - stats.min)*0.1;
        if margin <= 0 || isnan(margin), margin = 0.001; end
        xlim(ax, [stats.min - margin, stats.max + margin]);
    end

    title(ax, [plot_name ' (Zoomed)']);
    xlabel(ax, x_label);
    ylabel(ax, 'Count');
    grid(ax,'on');

    saveas(fig, fullfile(results_dir, ['hist' suffix '.svg']));
    fprintf('Saved hist%s.{svg}\n', suffix);
    fprintf(fileID,'Saved hist%s.{svg}\n', suffix);
    close(fig);
end

function plot_z_score_distributions(z_data, z_stats, results_dir, fileID)
    if isempty(z_data) || ~isfinite(z_stats.std)
        return;
    end

    z_mean = z_stats.mean;
    z_std  = z_stats.std;
    if z_std == 0, z_std = 1e-6; end
    SE   = z_stats.sem;
    df   = max(z_stats.N - 1, 1);
    tcrit= tinv(0.975, df);
    CI   = [z_mean - tcrit*SE, z_mean + tcrit*SE];

    min_val = z_mean - 4*z_std;
    max_val = z_mean + 4*z_std;
    if ~isfinite(min_val), min_val = -3; end
    if ~isfinite(max_val), max_val = 3;  end
    if min_val == max_val, min_val = min_val-1; max_val = max_val+1; end

    x   = linspace(min_val, max_val, 1000);
    pdf = normpdf(x, z_mean, z_std);

    fig = figure('Name','Z-Score Distribution','Visible','off', ...
                 'Position',[100,100,1000,600]);
    ax = gca; hold(ax,'on');
    plot(ax, x, pdf, 'LineWidth',2);
    yl = ylim(ax);
    line(ax, [z_mean z_mean], yl, 'LineStyle','--','LineWidth',1.5);
    line(ax, [CI(1) CI(1)], yl, 'LineWidth',2);
    line(ax, [CI(2) CI(2)], yl, 'LineWidth',2);
    ylim(ax, yl);
    title(ax, 'Temporal Reliability: Normal Approximation of Fisher z-transformed correlation values ');
    xlabel(ax, 'Fisher z-transformed Correlation Values'); ylabel(ax, 'Probality Density Function (PDF0');
    grid(ax,'on'); hold(ax,'off');

    saveas(fig, fullfile(results_dir, 'z_score_distributions.svg'));
    fprintf('Saved z_score_distributions.{svg}\n');
    fprintf(fileID,'Saved z_score_distributions.{svg}\n');
    close(fig);
end


