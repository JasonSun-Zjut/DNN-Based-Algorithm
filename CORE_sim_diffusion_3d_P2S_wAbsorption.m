
function [ nRx_wout_noise, n_destroy ] = CORE_sim_diffusion_3d_P2S_wAbsorption( ...
   tx_timeline, ...% Symbol sequence to modulate and transmit
   mol_type_cnt, ...% Molecule Type Count
   tx_node,    ...% Tx node properties
   rx_node,    ...% Rx node properties
   env_params, ...% Environment properties
   sim_params )   % Simulation parameters

p_react                 = rx_node.p_react;

if (p_react == 1)
   [ nRx_wout_noise, n_destroy ] = perfect_absorption(tx_timeline, mol_type_cnt, tx_node, rx_node, env_params, sim_params );
elseif (p_react == 2)
   [ nRx_wout_noise, n_destroy ] = imperfect_absorption(tx_timeline, mol_type_cnt, tx_node, rx_node, env_params, sim_params );
elseif (p_react == 3)
   [ nRx_wout_noise, n_destroy ] = mobile_perfect_absorption(tx_timeline, mol_type_cnt, tx_node, rx_node, env_params, sim_params );
else
    [ nRx_wout_noise, n_destroy ] = mobile_perfect_passive(tx_timeline, mol_type_cnt, tx_node, rx_node, env_params, sim_params );
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% different receiver model (passive, perfect_absorb, partical_absorb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ n_Rx, n_destroy ] = mobile_perfect_passive( ...
   tx_timeline, ...% Symbol sequence to modulate and transmit
   mol_type_cnt, ...
   tx_node,    ...% Tx node properties
   rx_node,    ...% Rx node properties
   env_params, ...% Environment properties
   sim_params )
rx_r_inMicroMeters      = rx_node.r_inMicroMeters;

tx_emission_point       = tx_node.emission_point;

D                       = env_params.D_inMicroMeterSqrPerSecond;
D_tx                    = tx_node.D_inMicroMeterSqrPerSecond;
D_rx                    = rx_node.D_inMicroMeterSqrPerSecond;
destruction_limit_sq    = env_params.destruction_limit^2;
delta_t                 = sim_params.delta_t;
alpha = 1.5;
% First find the number of simulation steps
sim_step_cnt = size(tx_timeline,2);

rx_membrane_sq = (rx_r_inMicroMeters)^2;


% Rx timeline Records the number of molecules release from Tx_1 at RECEIVER at each time step 
% Added DIFFERENT molecule TYPES by considering each row as another molecule type
n_Rx_wout_noise = zeros (mol_type_cnt, sim_step_cnt);

n_destroy = zeros (1, sim_step_cnt);
% mol released from tx_1
mol_position1 = zeros(0,3);
mol_type1 = ones(0,1);
for t=1:sim_step_cnt
   % Standard deviation of step size of molecules movement N(0,sigma)
    sigma_molecules = (2 * alpha * D*(t * delta_t)^(alpha-1) * delta_t)^0.5;
    % Standard deviation of step size of tx movement N(0,sigma)
    sigma_tx = (2 * alpha * D_tx*(t * delta_t)^(alpha-1) * delta_t)^0.5;
    % Standard deviation of step size of rx movement N(0,sigma)
    sigma_rx = (2 * alpha *D_rx*(t * delta_t)^(alpha-1) * delta_t)^0.5;
    
   % Check for Emission for EACH MOL_TYPE
   num_release = tx_timeline(:, t);
   if (sum(num_release) > 0)
      % Add new molecules to environment before moving them
      for ii=1:mol_type_cnt
          if (num_release(ii) > 0)
              mol_position1 = [ mol_position1 ; repmat(tx_emission_point(:,:), num_release(ii), 1) ];
              %mol_type1 = [ mol_type1; repmat([ii], num_release(ii), 1) ];
          end
      end
   end
    % Propagate the molecules via diffusion
    mol_displace = normrnd (0, sigma_molecules, size(mol_position1,1), 3);
    mol_displace_tx = normrnd(0, sigma_tx, size(tx_emission_point,1),3);
    mol_displace_rx = normrnd(0, sigma_rx, size(rx_node.center,1),3);
    tx_emission_point = tx_emission_point + mol_displace_tx;
    mol_position2 =  mol_position1 + mol_displace;
    rx_node.center = rx_node.center + mol_displace_rx;
    % Evaluate molecule from Tx_1 distance to Rx1 and Rx2
    dist_sq_2_rcv = sum(bsxfun(@minus, mol_position2, rx_node.center(1,:)).^2, 2);

    keep_mask = dist_sq_2_rcv < destruction_limit_sq;
    %the index of molecules need to be destroy%
    n_destroy(t) = n_destroy(t) + nnz(~keep_mask);
    %keep the ones indicated by the destruction mask (very far molecules are eliminated for efficiency)
    mol_position2 = mol_position2(keep_mask, :);
    dist_sq_2_rcv = dist_sq_2_rcv(keep_mask, :);
    %mol_type1 = mol_type1(keep_mask);
    %outside the membrane (continues its life)
    outside_membrane_mask_1 = dist_sq_2_rcv > rx_membrane_sq;
    %mol_type_mask = zeros(size(mol_type1,1), mol_type_cnt);
    %{
    for ii=1:mol_type_cnt
        mol_type_mask(:,ii) = (mol_type1 == ii);
    end
    %}
    %reception (hit)
    for ii=1:mol_type_cnt
        n_Rx_wout_noise(ii, t) = n_Rx_wout_noise(ii, t) + nnz(~outside_membrane_mask_1);
    end
    %absorption
    %{
    mol_position2 = mol_position2(outside_membrane_mask_1, :);
    mol_type1 = mol_type1(outside_membrane_mask_1);
    %}
    %passive
    mol_position1 = mol_position2;
end
    n_Rx = n_Rx_wout_noise;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% different receiver model (absorption, perfect_absorb, partical_absorb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ n_Rx_wout_noise, n_destroy ] = mobile_perfect_absorption( ...
   tx_timeline, ...% Symbol sequence to modulate and transmit
   mol_type_cnt, ...
   tx_node,    ...% Tx node properties
   rx_node,    ...% Rx node properties
   env_params, ...% Environment properties
   sim_params )
rx_r_inMicroMeters      = rx_node.r_inMicroMeters;

tx_emission_point       = tx_node.emission_point;

D                       = env_params.D_inMicroMeterSqrPerSecond;
D_tx                    = tx_node.D_inMicroMeterSqrPerSecond;
D_rx                    = rx_node.D_inMicroMeterSqrPerSecond;
destruction_limit_sq    = env_params.destruction_limit^2;
delta_t                 = sim_params.delta_t;

% First find the number of simulation steps
sim_step_cnt = size(tx_timeline,2);

% Standard deviation of step size of molecules movement N(0,sigma)
sigma_molecules = (2*D*delta_t)^0.5;
% Standard deviation of step size of tx movement N(0,sigma)
sigma_tx = (2*D_tx*delta_t)^0.5;
% Standard deviation of step size of rx movement N(0,sigma)
sigma_rx = (2*D_rx*delta_t)^0.5;

rx_membrane_sq = (rx_r_inMicroMeters)^2;


% Rx timeline Records the number of molecules release from Tx_1 at RECEIVER at each time step 
% Added DIFFERENT molecule TYPES by considering each row as another molecule type
n_Rx_1_wout_noise = zeros (mol_type_cnt, sim_step_cnt);
n_Rx_2_wout_noise = zeros (mol_type_cnt, sim_step_cnt);

n_destroy = zeros (1, sim_step_cnt);
% mol released from tx_1
mol_position1 = zeros(0,3);
mol_type1 = ones(0,1);
for t=1:sim_step_cnt
   % Check for Emission for EACH MOL_TYPE
   num_release = tx_timeline(:, t);
   if (sum(num_release) > 0)
      % Add new molecules to environment before moving them
      for ii=1:mol_type_cnt
          if (num_release(ii) > 0)
              mol_position1 = [ mol_position1 ; repmat(tx_emission_point(1,:), num_release(ii), 1) ];
              mol_type1 = [ mol_type1; repmat([ii], num_release(ii), 1) ];
          end
      end
   end
    % Propagate the molecules via diffusion
    mol_displace = normrnd (0, sigma_molecules, size(mol_position1,1), 3);
    mol_displace_tx = normrnd(0, sigma_tx, size(tx_emission_point,1),3);
    mol_displace_rx = normrnd(0, sigma_rx, size(rx_node.center,1),3);
    tx_emission_point = tx_emission_point + mol_displace_tx;
    mol_position2 =  mol_position1 + mol_displace;
    rx_node.center = rx_node.center + mol_displace_rx;
    % Evaluate molecule from Tx_1 distance to Rx1 and Rx2
    dist_sq_2_rcv = sum(bsxfun(@minus, mol_position2, rx_node.center(1,:)).^2, 2);
    dist_sq_2_rcv_2 = sum(bsxfun(@minus, mol_position2, rx_node.center(2,:)).^2, 2);
    % Evaluate molecule from Tx_2 distance to Rx1 and Rx2

    keep_mask = dist_sq_2_rcv < destruction_limit_sq & dist_sq_2_rcv_2 < destruction_limit_sq;
    %the index of molecules need to be destroy%
    n_destroy(t) = n_destroy(t) + nnz(~keep_mask);
    %keep the ones indicated by the destruction mask (very far molecules are eliminated for efficiency)
    mol_position2 = mol_position2(keep_mask, :);
    dist_sq_2_rcv = dist_sq_2_rcv(keep_mask, :);
    dist_sq_2_rcv_2 = dist_sq_2_rcv_2(keep_mask, :);
    mol_type1 = mol_type1(keep_mask);
    %outside the membrane (continues its life)
    outside_membrane_mask_1 = dist_sq_2_rcv > rx_membrane_sq;
    outside_membrane_mask_2 = dist_sq_2_rcv_2 > rx_membrane_sq;
    mol_type_mask = zeros(size(mol_type1,1), mol_type_cnt);
    for ii=1:mol_type_cnt
        mol_type_mask(:,ii) = (mol_type1 == ii);
    end
   
    %reception (hit)
    for ii=1:mol_type_cnt
        n_Rx_1_wout_noise(ii, t) = n_Rx_1_wout_noise(ii, t) + nnz(~outside_membrane_mask_1 & mol_type_mask(:,ii));
        n_Rx_2_wout_noise(ii, t) = n_Rx_2_wout_noise(ii, t) + nnz(~outside_membrane_mask_2 & mol_type_mask(:,ii));
    end
    %absorption
    mol_position2 = mol_position2(outside_membrane_mask_1 & outside_membrane_mask_2, :);
    mol_type1 = mol_type1(outside_membrane_mask_1 & outside_membrane_mask_2);
    %passive
    mol_position1 = mol_position2;
end
    n_Rx_wout_noise = [n_Rx_1_wout_noise, n_Rx_2_wout_noise];
end


function [ nRx_wout_noise, n_destroy ] = perfect_absorption( ...
   tx_timeline, ...% Symbol sequence to modulate and transmit
   mol_type_cnt, ...
   tx_node,    ...% Tx node properties
   rx_node,    ...% Rx node properties
   env_params, ...% Environment properties
   sim_params )   % Simulation parameters

rx_r_inMicroMeters      = rx_node.r_inMicroMeters;

tx_emission_point       = tx_node.emission_point;

D                       = env_params.D_inMicroMeterSqrPerSecond;
D_rx                    = rx_node.D_inMicroMeterSqrPerSecond;
destruction_limit_sq    = env_params.destruction_limit^2;

delta_t                 = sim_params.delta_t;


% First find the number of simulation steps
sim_step_cnt = size(tx_timeline,2);

% Standard deviation of step size of molecules movement N(0,sigma)
sigma_molecules = (2*D*delta_t)^0.5;
% Standard deviation of step size of rx movement N(0,sigma)
sigma_rx = (2*D_rx*delta_t)^0.5;

rx_membrane_sq = (rx_r_inMicroMeters)^2;


% Rx timeline Records the number of molecules at RECEIVER at each time step 
% Added DIFFERENT molecule TYPES by considering each row as another molecule type
nRx_wout_noise = zeros (mol_type_cnt, sim_step_cnt);

n_destroy = zeros (1, sim_step_cnt);

mol_position1 = zeros(0,3);
mol_type1 = ones(0,1);
for t=1:sim_step_cnt
   % Check for Emission for EACH MOL_TYPE
   num_release = tx_timeline(:, t);
   
   if (sum(num_release) > 0)
      % Add new molecules to environment before moving them
      for ii=1:mol_type_cnt
          if (num_release(ii) > 0)
              mol_position1 = [ mol_position1 ; repmat(tx_emission_point, num_release(ii), 1) ];
              mol_type1 = [mol_type1; repmat([ii], num_release(ii), 1)];
          end
      end
   end
   
    % Propagate the molecules via diffusion
    mol_displace = normrnd (0, sigma_molecules, size(mol_position1,1), 3);
    mol_displace_rx = normrnd(0,sigma_rx,size(rx_node.center,1),3);
    mol_position2 =  mol_position1 + mol_displace;
    rx_node.center = rx_node.center + mol_displace_rx;
    % Evaluate distance to Rx
    dist_sq_2_rcv = sum(bsxfun(@minus, mol_position2, rx_node.center).^2, 2);
    keep_mask = dist_sq_2_rcv < destruction_limit_sq;
    %the index of molecules need to be destroy%
    n_destroy(t) = n_destroy(t) + nnz(~keep_mask);
    %keep the ones indicated by the destruction mask (very far molecules are eliminated for efficiency)
    mol_position2 = mol_position2(keep_mask, :);
    dist_sq_2_rcv = dist_sq_2_rcv(keep_mask, :);
    mol_type1 = mol_type1(keep_mask);
    
    %outside the membrane (continues its life)
    outside_membrane_mask = dist_sq_2_rcv > rx_membrane_sq;
    
    mol_type_mask = zeros(size(mol_type1,1), mol_type_cnt);
    for ii=1:mol_type_cnt
        mol_type_mask(:,ii) = (mol_type1 == ii);
    end
    
    %reception (hit)
    for ii=1:mol_type_cnt
        nRx_wout_noise(ii, t) = nRx_wout_noise(ii, t) + nnz(~outside_membrane_mask & mol_type_mask(:,ii));
    end
    %nnz(~outside_membrane_mask & mol_type_mask(:,ii))
    %keep the ones indicated by the outside membrane mask
    mol_position2 = mol_position2(outside_membrane_mask, :);
    mol_type1 = mol_type1(outside_membrane_mask);
    mol_position1 = mol_position2;
end

% fprintf(1, '\n Destroyed (%d) Molecule(s) for Efficiency', sum(n_destroy));

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ nRx_wout_noise, n_destroy ] = imperfect_absorption( ...
   tx_timeline, ...% Symbol sequence to modulate and transmit
   mol_type_cnt, ...
   tx_node,    ...% Tx node properties
   rx_node,    ...% Rx node properties
   env_params, ...% Environment properties
   sim_params )   % Simulation parameters

p_react                 = rx_node.p_react;
rx_r_inMicroMeters      = rx_node.r_inMicroMeters;

tx_emission_point       = tx_node.emission_point;

D                       = env_params.D_inMicroMeterSqrPerSecond;
destruction_limit_sq    = env_params.destruction_limit^2;

delta_t                 = sim_params.delta_t;

% First find the number of simulation steps
sim_step_cnt = size(tx_timeline,2);

% Standard deviation of step size of movement N(0,sigma)
sigma = (2*D*delta_t)^0.5;


rx_membrane_sq = (rx_r_inMicroMeters)^2;

% Rx timeline Records the number of molecules at RECEIVER at each time step 
% Added DIFFERENT molecule TYPES by considering each row as another molecule type
nRx_wout_noise = zeros (mol_type_cnt, sim_step_cnt);

n_destroy = zeros (1, sim_step_cnt);

mol_position1 = zeros(0,3); 
mol_type1 = ones(0,1);
for t=1:sim_step_cnt
   % Check for Emission for EACH MOL_TYPE
   num_release = tx_timeline(:, t);
   
   if (sum(num_release) > 0)
      % Add new molecules to environment before moving them
      for ii=1:mol_type_cnt
          if (num_release(ii) > 0)
              mol_position1 = [ mol_position1 ; repmat(tx_emission_point, num_release(ii), 1) ];
              mol_type1 = [mol_type1; repmat([ii], num_release(ii), 1)];
          end
      end
   end
   
    % Propagate the molecules via diffusion
    mol_displace = normrnd (0, sigma, size(mol_position1,1), 3);
    mol_position2 =  mol_position1 + mol_displace;
    
    % Evaluate distance to Rx
    dist_sq_2_rcv = sum(bsxfun(@minus, mol_position2, rx_node.center).^2, 2);
    
    keep_mask = dist_sq_2_rcv < destruction_limit_sq;
    n_destroy(t) = n_destroy(t) + nnz(~keep_mask);
    %keep the ones indicated by the destruction mask (very far molecules are eliminated for efficiency)
    mol_position1 = mol_position1(keep_mask, :); % FOR ROLL BACK of ~preact MOLECULES
    mol_position2 = mol_position2(keep_mask, :);
    dist_sq_2_rcv = dist_sq_2_rcv(keep_mask, :);
    mol_type1 = mol_type1(keep_mask);
    
    %outside the membrane (continues its life)
    outside_membrane_mask = dist_sq_2_rcv > rx_membrane_sq;
    
    mol_type_mask = zeros(size(mol_type1,1), mol_type_cnt);
    for ii=1:mol_type_cnt
        mol_type_mask(:,ii) = (mol_type1 == ii);
    end
    
    %react random for all "TODO: Optimize this point"
    react_mask = random('Uniform', 0, 1, size(outside_membrane_mask)) < p_react;
    
    %reception (hit)
    for ii=1:mol_type_cnt
        nRx_wout_noise(ii, t) = nRx_wout_noise(ii, t) + nnz(~outside_membrane_mask & mol_type_mask(:,ii) & react_mask);
    end
    
    %keep the ones indicated by the outside membrane mask or nonReacting from mol_position1
    mol_position2 = [ mol_position2(outside_membrane_mask, :); mol_position1(~outside_membrane_mask & ~react_mask, :)];
    mol_type1 = [ mol_type1(outside_membrane_mask) ; mol_type1(~outside_membrane_mask & ~react_mask) ];

    mol_position1 = mol_position2;
end

end

