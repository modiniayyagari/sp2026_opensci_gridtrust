%% ENF Signature Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% This script implements the signal processing pipeline for the application 
% "Geolocation for Semiconductor Chips using Electric Network Frequency (ENF) Signatures".
% The primary goal is to validate the feasibility of extracting ENF signatures from 
% DC-powered hardware by comparing a sensed trace against a ground-truth
% reference.
%
% --- ANALYSIS OVERVIEW ---
% 1.  INPUTS: The script takes two time-synchronized inputs:
%     - A 'ground-truth reference trace' captured from the AC mains. (.wav file format)
%     - A 'sensed trace' captured from the experimental FPGA board.  (.wav file format)
%       This is the ambient EM trace.
%
% 2.  SPECTROGRAM GENERATION: The Short-Time Fourier Transform (STFT) is
%     applied to both traces to generate high-resolution spectrograms,
%     visualizing the harmonic content over time.
%
% 3.  ENF ESTIMATION: A weighted average PMF is used to extract the instantaneous ENF 
%     signature from each spectrogram for EM traces.
%
% 4.  OUTPUTS & CORRELATION: Finally, the script compares the ENF signature
%     from the sensed trace against the signature from the ground-truth 
%     reference. It calculates the Pearson correlation coefficient and 
%     generates the final temporally aligned plots to visually and quantitatively 
%     assess the match.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENF Analysis Configuration Flags:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input trace file information
% This script is designed to analyze ONE pair of traces from ONE specific
% directory (a "leaf folder") at a time. Manually update the variables below for EACH run:
%
% 1. `file_path`: Full relative path ending in a slash ("/") of the folder containing the trace pair.
% 2. `file_1_name`: File name of the AC mains reference trace without file extension.
% 3. `file_2_name`: File name of the sensed trace without file extension.
%
% --- EXAMPLES ---
%
% % Example 1: Analyzing a 'TREND' experiment trace (WK01/WEND/SUN/EVEN/SET2/T01)
% file_path = "../exp_inputs/TREND/WK01/WEND/SUN/EVEN/SET2/T01/";
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";
%
% % Example 2: Analyzing an 'FPGA' experiment trace (FPGA/SAKU/)
% file_path = "../exp_inputs/FPGA/SAKU/";
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";
%
% % Example 3: Analyzing a 'MULTI' (US_60) experiment trace (MULTI/US_60/AUG/WED/T02)
% % Note: File names are different in this folder.
% file_path = "../exp_inputs/MULTI/US_60/AUG/WED/T02/";
% file_1_name = "mains_pow_trace_ac_egrid_citya_lab";
% file_2_name = "fpga_em_trace_dc_egrid_citya_lab";
%
% % Example 4: Analyzing a 'MULTI' (DE_50) experiment trace (MULTI/DE_50/AUG/WED/T01)
% file_path = "../exp_inputs/MULTI/DE_50/AUG/WED/T01/";
% file_1_name = "mains_pow_trace_ac_citya_lab";
% file_2_name = "fpga_em_trace_dc_citya_lab";
%
% % Example 5: Analyzing an 'SRV_L' experiment trace (SRV_L/SPSU/DAY1/T05)
% file_path = "../exp_inputs/SRV_L/SPSU/DAY1/T05/";
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";

% Set the 3 variables below to point to the *single leaf folder* you want to analyze.
file_path = "../exp_inputs/FPGA/SAKU/"; 
file_1_name = "mains_pow_trace_ac";
file_2_name = "fpga_em_trace_dc";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrogram Computation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental power grid frequency settings
nominal_freq_arr      = [50 60];

% --- Multi-Location Experiment Note (DE_50 / US_60) ---
% When analyzing traces from the 'MULTI' directory, you MUST manually
% update the `nominal_freq_*` variables below to match the fundamental
% frequency of the grid from where the trace was recorded.
% - For US_60 (United States) traces: set nominal_freq = nominal_freq_arr(2); (60 Hz)
% - For DE_50 (Germany) traces:     set nominal_freq = nominal_freq_arr(1); (50 Hz)

%For reference trace
nominal_freq_1        = nominal_freq_arr(2);        
harmonics_arr_1       = [1:7]*nominal_freq_1;

%For sensed trace 
nominal_freq_2        = nominal_freq_arr(2);        
harmonics_arr_2       = [1:7]*nominal_freq_2;

% STFT compute param settings
frame_size_arr      = [1:12]*1000;
frame_size          = frame_size_arr(8);                %8000ms window
nfft_arr            = 2.^[10:20];
nfft                = nfft_arr(6);                      %2^15 = 32768 pts
overlap_size_arr    = 0:0.1:0.9;
overlap_size        = overlap_size_arr(1)*frame_size;   %non-overlapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Estimation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency to be estimated
trace_1_est_freq = harmonics_arr_1(1);                    %1st harmonic= 60Hz    
trace_2_est_freq = harmonics_arr_2(1);                    %1st harmonic= 60Hz

%Frequency estimation method options: 
%1: weighted average (pmf power = 3)
%2: spectrum combining (quad interp)
trace_1_freq_est_method = 1;   %For reference trace                         
trace_1_freq_est_spec_comb_harmonics = [60 120 180 240 300 360 420];

trace_2_freq_est_method = 1;   %For sensed trace                        
trace_2_freq_est_spec_comb_harmonics = [60 120 180 240 300 360 420];

tempo_align_corr_vals = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Plotting Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The colormap for the figures is set to jet.
set(0,'DefaultFigureColormap', jet)

%Titles for MATLAB plots
trace_1_plot_title = "Reference AC Mains Power Trace";
trace_2_plot_title = "Ambient EM Trace";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define .wav file paths for the rest of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_file_1_path_wav = file_path + file_1_name + ".wav";
full_file_2_path_wav = file_path + file_2_name + ".wav";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENF Analysis for the complete trace (calls the pre-compiled executable which contains proprietory code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handles spectrogram computation, ENF extraction, ENF
% matching and plotting all the figures

fprintf('\nINFO: Performning ENF Analysis on reference and sensed traces ...\n');
matched_corr_val_full = proc_enf_analysis(full_file_1_path_wav, full_file_2_path_wav, ...
                                          nfft, frame_size, overlap_size, ...
                                          harmonics_arr_1, nominal_freq_1, harmonics_arr_2, nominal_freq_2,...
                                          trace_1_freq_est_method, trace_1_est_freq, trace_1_freq_est_spec_comb_harmonics, ...
                                          trace_2_freq_est_method, trace_2_est_freq, trace_2_freq_est_spec_comb_harmonics, ...
                                          trace_1_plot_title, trace_2_plot_title, true);
                                      
tempo_align_corr_vals{end+1} = matched_corr_val_full;
fprintf('\nINFO: Pearsons Correlation Coefficient for the temporally aligned reference and sesned traces (complete trace): %.8f \n', cell2mat(tempo_align_corr_vals(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

