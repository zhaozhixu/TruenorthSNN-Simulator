% simulator for Truenorth SNN architecture
% top block
uint_bit_width = 8;
int_bit_width = 9;
PRNceil = 2 ^ uint_bit_width - 1;
max_rate = 1000;

fprintf('Truenorth SNN Architecture Simulator V0.02\n');
fprintf('[1] Configure a new network\n');
fprintf('[2] Use an old configuration file(.mat)\n');
option = input('');
if option == 1
  fprintf('Axon number :\n');
  axon_num = input('');
  fprintf('Neuron number :\n');
  neu_num = input('');
  fprintf('Time interval :\n');
  time_interval = input('');

  V = zeros(1, neu_num); %membrane potential
  w = zeros(axon_num, neu_num); %synapse connection
  s = zeros(axon_num, neu_num); %synaptic weight
  b = zeros(axon_num, neu_num); %synapse mode
  c = zeros(1, neu_num); %leak mode
  R = zeros(1, neu_num); %reset voltage
  lamda = zeros(1, neu_num); %leak weight
  alpha = zeros(1, neu_num); %positive threshold
  beta = zeros(1, neu_num); %negative threshold
  gamma = zeros(1, neu_num); %reset mode
  kapa = zeros(1, neu_num); %negative thresh: reset or saturate
  epsilon = zeros(1, neu_num); %leak-reversal flag
  M = zeros(1, neu_num); %threshold PRN mask
  rho_s = zeros(axon_num, neu_num); %synaptic PRN
  rho_l = zeros(1, neu_num); %leak PRN
  rho_t = zeros(1, neu_num); %threshold PRN

  for i = 1:neu_num
    fprintf('Synapse connection vector of neuron %d with %d axons:\n', i, axon_num);
    w(:, i) = input('');
    snp_num = sum(w(:, i));
    index = find(w(:, i) == 1);
    fprintf('Synapse mode vector of neuron %d''s %d synapses:\n', i, snp_num);
    b(index, i) = input('');
    fprintf('Synapse weight vector of neuron %d''s %d synapses:\n', i, snp_num);
    s(index, i) = input('');
  end
  fprintf('Leak mode vector of %d neurons:\n', neu_num);
  c = input('');
  fprintf('Leak weight vector of %d neurons:\n', neu_num);
  lamda = input('');
  fprintf('Leak-reversal flag vector of %d neurons:\n', neu_num);
  epsilon = input('');
  fprintf('Reset mode vector of %d neurons:\n', neu_num);
  gamma = input('');
  fprintf('Reset voltage vector of %d neurons:\n', neu_num);
  R = input('');
  fprintf('Positive threshold vector of %d neurons:\n', neu_num);
  alpha = input('');
  fprintf('Negative threshold vector of %d neurons:\n', neu_num);
  beta = input('');
  fprintf('Threshold PRN mask vector of %d neurons(0 to %d):\n', neu_num, PRNceil);
  M = input('');
  fprintf('Negative thresh mode vector of %d neurons:\n', neu_num);
  kapa = input('');

elseif option == 2
  fprintf('Input the configuration file name:\n');
  path = input('', 's');
  load(path);
else
  error('Please input a valid option!\n');
end

axon = zeros(axon_num, time_interval); %axon input
spike = zeros(neu_num, time_interval); %spike output
in = zeros(axon_num, max_rate); %input of spike rate
out = zeros(neu_num, max_rate); %output of spike rate
fprintf('Input the coefficient vector of spike rate expression for every axon.\n');
fprintf('Linear expressions only. For example, if the expression is x+2, input [1 2]. \n');
for i = 1:axon_num
  fprintf('Spike rate expression coefficient vector of axon %d:\n', i);
  expr = input('');
  in(i, :) = round(expr(1) * (1:max_rate) + expr(2));
end
fprintf('Computing...');

for j = 1:max_rate
  axon = rate_input(in(:, j), time_interval);
  for i = 1:time_interval
    index = find(b == 1);
    if ~isempty(index)
      rho_s(index) = round(rand(1, length(index)) * PRNceil);
    end
    index = find(c == 1);
    if(~isempty(index))
      rho_l(index) = round(rand(1, length(index)) * PRNceil);
    end
    index = find(M ~= 0);
    if ~isempty(index)
      rho_t(index) = round(rand(1, length(index)) * PRNceil);
    end
    V = V + sum(expand(axon(:, i), neu_num).* w.* ((1 - b).* s +...
                                                   b.* bin_cmp(s, rho_s).* sign(s)));

    V = V + ((1 - epsilon + epsilon.* sign(V)).* ((1 - c).* lamda) +...
             c.* bin_cmp(lamda, rho_l).* sign(lamda));
    eta = bitand(rho_t, M); %masked threshold PRN
    spike(:, i) = transpose(V >= alpha + eta);
    index = find(spike(:, i) == 1);
    if ~isempty(index)
      V(index) = delta_fun(gamma(index)).* R(index) +...
                 delta_fun(gamma(index) - 1).* (V(index) - (alpha(index) + eta(index))) +...
                 delta_fun(gamma(index) - 2).* V(index);
    end
    index = find(V <= -(beta.* kapa + (beta + eta).* (1 - kapa)) == 1);
    if ~isempty(index)
      V(index) = -beta(index).* kapa(index) +...
                 (-delta_fun(gamma(index)).* R(index) +...
                  delta_fun(gamma(index) - 1).* (V(index) + (beta(index) + eta(index))) +...
                  delta_fun(gamma(index) - 2).* V(index)).* (1 - kapa(index));
    end
  end
  out(:, j) = sum(spike, 2) / (time_interval / max_rate);
end
plot(1:max_rate, in, 'g.', 1:max_rate, out, 'r.');
%legend('input', 'input', 'output', 'location', 'southoutside');
fprintf('finished.\n');
option2 = input('Save current configuration file?[y/n]', 's');
if option2 == 'y' || option2 == 'Y'
  save_path = input('Save path:', 's');
  save(save_path);
elseif option2 == 'n' || option2 == 'N'
else
  error('Please answer y or n.\n');
end
