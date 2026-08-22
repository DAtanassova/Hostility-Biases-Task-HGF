clearvars;

%% Specify paths

% Project Location 
% specify where the data for the model fitting is
% format: for each participant, a file in the format ("ppn_valence_input_pu.txt")
% should exist, where the first column is the valence of the stimulus
% (here: 1=hostile, 0=non-hostile) and the second is the perceptual input
% including ambiguity (0.95/0.05/0.05).
% additionally, there should be a file in the format ("ppn_response.txt")
% containing the choice on each trial (0/1) and the trialwise RTs (in ms)

project_loc = ['~/Library/CloudStorage/OneDrive-RadboudUniversiteit/Work Documents/Projects/Hostility Task Danique/']
input_data_location = [project_loc 'data/hgf model files'];
save_results_location = [project_loc 'hgf model results'];

%% Specify the subjectIDs

vp = {'1','2','3',...};

%% Load model space
% original tapas directory: download this locally and specify the path to
% where the folder is located
addpath '~/HGF Toolbox/HGF'
% custom models: download the custom models for the HB task and specify the
% path to where they're located
addpath '~/HGF Toolbox/HB models'

%% Specify the prior parameter configuration of the perceptual model
% here, we use the adapted tapas_ehgf_ar1_binary_pu model with its default
% parameter space
perc_model2_config = tapas_ehgf_ar1_binary_pu_config();
% you can change the prior parameter configuration in the following way:
% perc_model2_config.ommu(1,2) = -4; % changes the prior mean of om2=-4
% the same can be done with the response model (not shown here)

%% Run the model and simulate trialwise RTs
% specify how many simulations to run (n=100 recommended)
n_sim = 100;

% specify the structure where the simulated data will be kept
sim_y_mean_all = table();

for i = 1:length(vp)
    
    disp(strcat("model fitting: subj  ", num2str(i)));

    % load the input file for that participant
    u_file_path = fullfile(input_data_location, vp(i) + "_valence_input_pu.txt");  
    u = load(u_file_path);
    % load the choice+RT file for that participant
    file_path = fullfile(input_data_location, vp(i) + "_response.txt");  
    y = load(file_path);

    % fit the model: y=choices, u=input, perc_model2_config chosen
    % perceptual model; "tapas_logrt_linear_binary_val_config" chosen
    % response model; "tapas_quasinewton_optim_config" optimization
    % algorithm
    est1 = tapas_fitModel(y,...
                        u,...
                       perc_model2_config,... 
                       'tapas_logrt_linear_binary_val_config',...
                       'tapas_quasinewton_optim_config'); 


 
    all_fits_model1(i) = est1; 
    % save the structure (comprising the estimated parameters and
    % trialwise computational quantities for each participant)
    model1_file_path = fullfile(save_results_location, 'logrt_val_binary_ehgf_ar1_pu.mat');
    save(model1_file_path, 'all_fits_model1')
    est = est1;

    % Begin simulations
    n = n_sim;
    c = 1;

       while c <= n 
        % for each participant, simulate data with the chosen perceptual
        % model and the chosen response model. The individual parameters
        % are added through est.p_prc.p/est.p_obs.p.
        sim_temp = tapas_simModel(u,...
                         'tapas_ehgf_ar1_binary_pu',...
                         [est.p_prc.p],...
                         'tapas_logrt_linear_binary_val',...
                         [est.p_obs.p]);
                     
       sim_temp.vp = vp(i);
       sim(c) = sim_temp;
       
       c = c+1;
            
       end 

       % save all simulations for that subject in a struct
       simulated_data_struct{i} = sim;
                    
        for n = 1:n_sim
               sim_y(:,n) = sim(n).y; 

        end
        
        num_rows = size(y, 1);
        sim_y_mean = NaN(num_rows, 3);
          
        sim_y_mean(:,1) = mean(sim_y,2);
        sim_y_mean(:,2) = log(y(:,2));
        sim_y_mean(:,3) = 1:num_rows;
        sim_y_mean(:,4) = u(:,1); 
        sim_y_mean(:,5) = sim(7).y; 
        sim_y_mean(:,6) = u(:,2); 

        sim_y_mean_table = array2table(sim_y_mean);
        sim_y_mean_table(:, 7) = repmat(vp(i), num_rows, 1);
        sim_y_mean_table.Properties.VariableNames = {'mean_sim_RT', ...
            'og_RT', 'trial', 'target', 'random_sim_RT', 'val', 'subjID'};
        sim_y_mean_all = [sim_y_mean_all; sim_y_mean_table];
             
        % save the simulated data for all participants in a .mat file and a
        % separate .csv file
        sims_struct_file_path = fullfile(save_results_location, 'simulations_logrt_val_binary_ehgf_ar1_pu.mat');
        save(sims_struct_file_path, 'simulated_data_struct')
        sims_table_file_path = fullfile(save_results_location, 'simulations_logrt_val_binary_ehgf_ar1_pu.csv');
        writetable(sim_y_mean_all, sims_table_file_path, 'WriteVariableNames', true); 

end     

% When simulating data, you can check whether the original behavior trends
% (for exampke, faster RTs to hostile faces) are also preserved in the
% simulated data. 

%% Recover parameters
% Recover parameters to determine their identifiability. Parameters that
% are poorly recovered (e.g. correlation between true and recovered
% parameter is <0.5) should be excluded from analyses as they are not
% reliably identifiable.

% Load all simulations (generated in the step above)
all_sims = load([sims_struct_file_path]);

% Speficy the number of simulations that should be recovered (here: 10).
% The simulations will be chosen at random.
n_recovered = 10;

for i = 1:length(vp)
    
    disp(strcat("model recovery: subj  ", num2str(i)));

    u_file_path = fullfile(input_data_location, vp(i) + "_valence_input_pu.txt");  
    u = load(u_file_path);
    % load the choice+RT file for that participant
    file_path = fullfile(input_data_location, vp(i) + "_response.txt");  
    y = load(file_path);

    current_subj_data = all_sims.simulated_data_struct{1, i};
    num_elements = numel(current_subj_data);
    random_indices = randperm(num_elements, n_recovered);
    random_sims = current_subj_data(random_indices);

n = n_recovered;
c = 1;

       while c <= n 
           
        sim_y = exp(random_sims(c).y);
      % logrt model expects RT in the second column so double them
        sim_y(:,2) = sim_y;
           
        % refit the same combination of model to the simulated data (here:
        % simulated responses are stored in sim_y);
        rec_est_temp = tapas_fitModel(sim_y,...
                        u,...
                       perc_model2_config,... 
                       'tapas_logrt_linear_binary_val_config',...
                       'tapas_quasinewton_optim_config'); 
                   
        rec_est_temp.vp = vp(i);
        rec_est_temp.n_fits = n_recovered;
        
        rec_est(c) = rec_est_temp;

       c = c+1;
            
       end 
    
       % save all recovered fits for that subject in a struct and save it
       recovered_data_struct{i} = rec_est;
       rec_est_struct_file_path = fullfile(save_results_location, 'recovered_fits_logrt_val_binary_ehgf_ar1_pu.mat');
       save(rec_est_struct_file_path, 'recovered_data_struct')

end 


%% Extract the true parameters

% specify the name of the file to be saved
parameters_file_path = fullfile(save_parameters_location, 'estimated_parameters_logrt_val.csv');

% load all the recovered data
all_fits = load([model1_file_path]);

for i = 1:length(vp)
    
    current_subj_fits = all_fits.all_fits_model1(i);

    est_parameters(i,1) = current_subj_fits.p_prc.om(1,2); % om2: tonic volatility
    est_parameters(i,2) = current_subj_fits.p_prc.phi(1,2); % phi2: belief resetting
    est_parameters(i,3) = current_subj_fits.p_prc.al(1,1); % alpha: perceptual uncertainty
    est_parameters(i,4) = current_subj_fits.p_obs.be1(1,1); % condition effect
    est_parameters(i,5) = current_subj_fits.p_obs.be2(1,1); % surprise effect  
    est_parameters(i,6) = current_subj_fits.p_obs.be3(1,1); % interaction  
    est_parameters(i,7) = current_subj_fits.p_obs.ze(1,1); % RT noise  

    est_parameters_table = array2table(est_parameters);
    est_parameters(i,8) = str2double(vp{i});


end 

est_parameters_table.Properties.VariableNames = {'om2', 'phi2', 'alpha', ...
    'be1_val', 'be2_surp', 'be3_val_surp_int', 'ze', 'subjID'};

writetable(est_parameters_table, parameters_file_path, 'WriteVariableNames', true); 

%% Extract trialwise trajectories
% specify the name of the file to be saved
trajectories_file_path = fullfile(save_parameters_location, 'estimated_trajectories_logrt_val.csv');

% load all the recovered data
all_fits = load([model1_file_path]);

trajectories_all = table();

for i = 1:length(vp)
    
    current_subj_fits = all_fits.all_fits_model1(i);

    num_rows = size(current_subj_fits.y, 1);

    trajectories_subj(:,1) = current_subj_fits.traj.muhat(:,1);
    trajectories_subj(:,2) = current_subj_fits.traj.mu(:,2);
    trajectories_subj(:,3) = current_subj_fits.traj.mu(:,3);
    trajectories_subj(:,4) = current_subj_fits.traj.sahat(:,1);
    trajectories_subj(:,5) = current_subj_fits.traj.sa(:,2);
    trajectories_subj(:,6) = current_subj_fits.traj.sa(:,3);
    trajectories_subj(:,7) = current_subj_fits.traj.epsi(:,2);
    trajectories_subj(:,8) = current_subj_fits.traj.epsi(:,3);
    trajectories_subj(:,9) = current_subj_fits.u(:,1); % target
    trajectories_subj(:,10) = current_subj_fits.u(:,2); % valence
    trajectories_subj(:,11) = current_subj_fits.y(:,1); % choice
    trajectories_subj(:,12) = current_subj_fits.y(:,2); % RT
    trajectories_subj(:,13) = 1:num_rows;

    trajectories_table = array2table(trajectories_subj);
    trajectories_table(:, 14) = repmat(vp(i), num_rows, 1);
    trajectories_table.Properties.VariableNames = {'mu1hat', ...
            'mu2', 'mu3', 'sa1hat', 'sa2', 'sa3', 'epsi1',...
            'epsi2','target','valence','choice','RT','trial','subjID'};
    trajectories_all = [trajectories_all; trajectories_table];

end 

    
writetable(trajectories_all, trajectories_file_path, 'WriteVariableNames', true); 

%% Extract parameters: recovered and true

% specify all structures
model1_file_path = fullfile(save_results_location, 'logrt_val_binary_ehgf_ar1_pu.mat');
rec_est_struct_file_path = fullfile(save_results_location, 'recovered_fits_logrt_val_binary_ehgf_ar1_pu.mat');
parameters_file_path = fullfile(save_parameters_location, 'estimated_and_recovered_parameters_logrt_val.csv');

% load the true and recovered data
all_fits = load([model1_file_path]);
all_recs = load([rec_est_struct_file_path]);

for i = 1:length(vp)

    current_subj_fits = all_fits.all_fits_model1(i);
    current_subj_recs = all_recs.recovered_data_struct(i);

    for n = 1:n_recovered
              om2_sim(1,n) = current_subj_recs{1,1}(n).p_prc.om(1,2);
              phi2_sim(1,n) = current_subj_recs{1,1}(n).p_prc.phi(1,2);
              alpha_sim(1,n) = current_subj_recs{1,1}(n).p_prc.al(1,1);
              be1_sim(1,n) = current_subj_recs{1,1}(n).p_obs.be1(1,1);
              be2_sim(1,n) = current_subj_recs{1,1}(n).p_obs.be2(1,1); % surprise effect  
              be3_sim(1,n) = current_subj_recs{1,1}(n).p_obs.be3(1,1); % interaction  
              ze_sim(1,n) = current_subj_recs{1,1}(n).p_obs.ze(1,1); % noise  
        end
        
        om2_sim_mean = mean(om2_sim);
        phi2_sim_mean = mean(phi2_sim);
        alpha_sim_mean = mean(alpha_sim);
        be1_sim_mean = mean(be1_sim);
        be2_sim_mean = mean(be2_sim);
        be3_sim_mean = mean(be3_sim);
        ze_sim_mean = mean(ze_sim);

    rec_og_parameters(i,1) = current_subj_fits.p_prc.om(1,2);
    rec_og_parameters(i,2) = current_subj_fits.p_prc.phi(1,2);
    rec_og_parameters(i,3) = current_subj_fits.p_prc.al(1,1);
    rec_og_parameters(i,4) = current_subj_fits.p_obs.be1(1,1);
    rec_og_parameters(i,5) = current_subj_fits.p_obs.be2(1,1);
    rec_og_parameters(i,6) = current_subj_fits.p_obs.be3(1,1);
    rec_og_parameters(i,7) = current_subj_fits.p_obs.ze(1,1);
    rec_og_parameters(i,8) = om2_sim_mean;
    rec_og_parameters(i,9) = phi2_sim_mean;
    rec_og_parameters(i,10) = alpha_sim_mean;
    rec_og_parameters(i,11) = be1_sim_mean;
    rec_og_parameters(i,12) = be2_sim_mean;
    rec_og_parameters(i,13) = be3_sim_mean;
    rec_og_parameters(i,14) = ze_sim_mean;
    rec_og_parameters(i, 15) = str2double(vp{i});

    rec_og_parameters_table = array2table(rec_og_parameters);


end 

rec_og_parameters_table.Properties.VariableNames = {'om2_og','phi2_og',...
    'alpha_og','be1_og', 'be2_og', 'be3_og', 'zeta_og',...
    'om2_sim','phi2_sim', 'alpha_sim',...
    'be1_sim', 'be2_sim', 'be3_sim', 'zeta_sim', 'ppn'};
% save the data
writetable(rec_og_parameters_table, parameters_file_path, 'WriteVariableNames', true); 
