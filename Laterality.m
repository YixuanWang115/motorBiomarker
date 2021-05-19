% subName = {''};
p_beta_right = 0; p_beta_left = 0;
p_alpha_right = 0; p_alpha_left = 0;
laterality1 = zeros(20,1);laterality2 = zeros(20,1);
Fs = 128;
%%
for subInd = 1 : 20
    %% Extract Blocks
    signal = cell2mat(struct2cell(load(append(subName{subInd},'.mat'))));
    eventInf = xlsread(append(subName{subInd},'.csv'));
    eventTyp = eventInf(:,4);
    eventTime = round(Fs.*eventInf(:,1)); % Convert Time to index
    
    % Search based on event types (may need to change values)
    % Search for statring points of the blocks
    % Change codes based on left/right hand tasks
    % eventLoc = round(find(eventTyp==1));
    % eventLoc = find(eventTyp==2);
    rstLoc = round(find(eventTyp==3));
    task_all = {};
    for i = 1 : length(Ts)
        task_all_l{i} = signal(11,eventTime(eventLoc(i))+0.5*Fs:eventTime(eventLoc(i)) + 3.5*Fs);
        task_all_r{i} = signal(32,eventTime(eventLoc(i))+0.5*Fs:eventTime(eventLoc(i)) + 3.5*Fs);
    end
    rst_all ={};
    for i = 1 : length(Te)
        rst_all_l{i} = signal(11,eventTime(rstLoc(i))+0.5*Fs:eventTime(rstLoc(i))+3.5*Fs);
        rst_all_r{i} = signal(32,eventTime(rstLoc(i))+0.5*Fs:eventTime(rstLoc(i))+3.5*Fs);
    end 
    p_ERD_right = getERD(task_all_r,rst_all_r);
    p_ERD_left = getERD(task_all_l,rst_all_r);
    p_ERD_right = p_ERD_right-min([p_ERD_right,p_ERD_left]);
    p_ERD_left = p_ERD_left-min([p_ERD_right,p_ERD_left]);
    % Left hand
%          laterality1(subInd) = (mean(p_ERD_right)-mean(p_ERD_left))/(mean(p_ERD_right)+mean(p_ERD_left));
    % Right hand
%   laterality2(subInd) = (mean(p_ERD_left)-mean(p_ERD_right))/(mean(p_ERD_right)+mean(p_ERD_left));
end