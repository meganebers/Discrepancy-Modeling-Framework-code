function [idx_heel,idx_knee] = cycleID(t,x,steps)

% Identify jumps in data (heelstrike collisions)
stride = x(:,1);
stride = diff(stride);% differentiate to find heelstrike (discontinuity in states)

% figure,
% plot(stride,'Linewidth',[2])
% title('Heelstrike'),
% xlabel('Time'),ylabel('')

thresh = 0.5*max(stride);

stride(stride < thresh) = 0;

% plot(stride,'Linewidth',[2])
% title('Residuals')

idx_heel = find(stride ~= 0);
% idx_startend = [1;idx;length(x)];
% 
% numStrides = numel(idx)+1;

knee_condition = find(abs(x(:,2) - x(:,3)) == 0);
knee_collision = diff(knee_condition);
knee_collision(knee_collision < thresh) = 0;
idx_knee = find(knee_collision ~= 0);

check_heel = numel(idx_heel)+1;
check_knee = numel(idx_knee)+1;

if check_heel && check_knee == steps
    ;
else
    disp('Number of strides was not correctly identified OR data is too noisy');
    return;
end

% %extracts gait cycles for each stride without discontinuity
% x_stepExtract = [];
% %dx_stepExtract = [];
% t_stepExtract = [];
% for j = 1:length(idx_startend)-1
%     x_stepExtract{j} = x(idx_startend(j)+1:idx_startend(j+1),:);
%     %dx_stepExtract{j} = dx(idx_startend(j)+1:idx_startend(j+1),:);
%     t_stepExtract{j} = t(idx_startend(j)+1:idx_startend(j+1),:);
% end
% 
% %         x_temp = [];
% %         for k = 1:length(idx)
% %             x_temp{k} = x_stepExtract{k};
% %         end
% 
% t_s = cell2mat(t_stepExtract.');
% %dx_s = cell2mat(dx_stepExtract.');
% x_s = cell2mat(x_stepExtract.');