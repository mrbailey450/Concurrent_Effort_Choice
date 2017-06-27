function [DATA] = CEC_Bout_Pause_Dat(Var,TS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% - Var = rawlist.data
% - TS = trial_structure from EFFORT_DATA


% - - This is a function which will determine the duration of every lever
% press, the inter-response-times (IRT's) or Pauses, as well as the
% structure of the bouts of responding: Number of Bouts, Lenth of each Bout
% Ect


% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
n = length(Var);

%LLeverOn = 0027 \***
%RLeverOn= 0028  \***

%ForcedTrial = 3300
%ChoiceTrial = 3301

%PelletTrial = 3320
%SucroseTrial = 3321

%Pelletchoice = 0100\code for hold choice
%Sucrosechoice = 0101\code for press choice
%ForcedOpOut = 03330\opting out of forced trial

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


% - - - All just shit to allow me to know dimensions of TS - - - 
%=======================================================%
%tells how many start trial/end trials are there%
for i = 1:n
    num_start{1,i} = length(find(Var{1,i}(:,2)==111));
    num_end{1,i} = length(find(Var{1,i}(:,2)==112));
end

%Defines rows in Ru of Start, end, AND End Session%
for i = 1:n
    start_rows{1,i} = find(Var{1,i}(:,2)==111);
    end_rows{1,i} = find(Var{1,i}(:,2)==112);
    end_all{1,i} = find(Var{1,i}(:,2)==118);
end

%Abbrevating rows with start trial%
sr = start_rows;


%Subtracting 1 from start trial b/c sr-1 is end of trial%
for i = 1:n
    sr_m1{1,i} = sr{1,i}(:,1) - 1;
end

%removing first row%
for i = 1:n
    sr_m1{1,i}(1,:)=[];
end


%Adding position of session end to end of sr_m1, b/c thats the
%row # for the end of the last trial%
for i = 1:n
    sr_m1{1,i}(end+1,:) = end_all{1,i};
end

%Variable defining # of start and end trials%
for i = 1:n
    Sess_Info_num_trials{1,i} = [num_start{1,i};num_end{1,i}];
end
SI_nt = Sess_Info_num_trials;   %easier var name to work with%

%Easy name for how many trials are there%
for i = 1:n
    A{1,i} = SI_nt{1,i}(1,1);
end

%Making a TRIALS data structure: each row contains all time stamps
%and event codes from that trial number%

%NOTE: switched "i" here and now "j" = # subs (i.e.13)
for j = 1:n
    for i = 1:A{1,j}
        trial_structure{i,j} = Var{1,j}(sr{1,j}(i,1):sr_m1{1,j}(i,1),:);
    end
end
% - - - All just shit to allow me to know dimensions of TS - - -

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

TS_Lever_Zero = TS;

%Easier name to use
TSL = TS_Lever_Zero;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -
%Seting the start time of each trial to zero

            % Two Steps: A. and B. 

% A. 
%Figuring out the start time of each trial (Globally)
for i = 1:n
    for j = 1:A{1,i}
        TSL_Trial_Start_Global_Time(j,i) =  TSL{j,i}(1,1);
    end
end 

% B.
%Subtracting this time from all of the other times in the TSL struct within
%each trial - Zeroing everything so each trial begins at time zero

for i = 1:n
    for j = 1:A{1,i}
        TSL{j,i}(:,1) = TSL{j,i}(:,1) - TSL_Trial_Start_Global_Time(j,i);
    end
end 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -


% - Determining the row number that the lever comes out for each trial
for i = 1:n
    for j = 1:A{1,i}
        if length(find(TSL{j,i}(:,2) == 27)) == 0
            TSL_LeverOut_Rows{j,i} = find(TSL{j,i}(:,2) == 28);
        else
            TSL_LeverOut_Rows{j,i} = find(TSL{j,i}(:,2) == 27);
        end
    end
end

% - Determining the time that the lever comes out for each trial

for i = 1:n
    for j = 1:A{1,i}
        TSL_LeverOut_Time(j,i) = TSL{j,i}(TSL_LeverOut_Rows{j,i},1);
    end
end


% - -  Getting just the forced trials
for i = 1:n
    for j = 1:10
        TSL_Forced{j,i} = TSL{j,i};
    end
end

TSLF = TSL_Forced;

% - Determining the Press Rows
for i = 1:n
    for j = 1:10
        if length(find(TSL{j,i}(:,2) == 28)) == 0
            TSL_Any_Press_Rows{j,i} = find(TSLF{j,i}(:,2) == 1015);
        else
            TSL_Any_Press_Rows{j,i} = find(TSLF{j,i}(:,2) == 1016);
        end
    end
end



for i = 1:n
    P = 1;
    H = 1;
    for j = 1:10
        if length(find(TSLF{j,i}(:,2) == 3320)) == 0
            PRESS_STRUCT{P,i} = TSLF{j,i};
            P = P + 1;
        else
            HOLD_STRUCT{H,i} = TSLF{j,i};
            H = H + 1;
        end
    end
end

% - - Determining the Pellet Rows and Hold Rows

for i = 1:n
    for j = 1:5
        if length(find(PRESS_STRUCT{j,i}(:,2) == 28)) == 0
            Press_Down_Rows{j,i} = find(PRESS_STRUCT{j,i}(:,2) == 1015);
            Press_Up_Rows{j,i} = find(PRESS_STRUCT{j,i}(:,2) == 1017);
            Forced_Press_Lever_Out_Rows(j,i) = find(PRESS_STRUCT{j,i}(:,2) == 27);
        else
            Press_Down_Rows{j,i} = find(PRESS_STRUCT{j,i}(:,2) == 1016);
            Press_Up_Rows{j,i} = find(PRESS_STRUCT{j,i}(:,2) == 1018);
            Forced_Press_Lever_Out_Rows(j,i) = find(PRESS_STRUCT{j,i}(:,2) == 28);
        end
    end
end


for i = 1:n
    for j = 1:5
        if length(find(HOLD_STRUCT{j,i}(:,2) == 28)) == 0
            Hold_Down_Rows{j,i} = find(HOLD_STRUCT{j,i}(:,2) == 1015);
            Hold_Up_Rows{j,i} = find(HOLD_STRUCT{j,i}(:,2) == 1017);
            Forced_Hold_Lever_Out_Rows(j,i) = find(HOLD_STRUCT{j,i}(:,2) == 27);
        else
            Hold_Down_Rows{j,i} = find(HOLD_STRUCT{j,i}(:,2) == 1016);
            Hold_Up_Rows{j,i} = find(HOLD_STRUCT{j,i}(:,2) == 1018);
            Forced_Hold_Lever_Out_Rows(j,i) = find(HOLD_STRUCT{j,i}(:,2) == 28);
        end
    end
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for i = 1:n
    for j = 1:5
        PRESS_DOWN_TIMES{j,i} = PRESS_STRUCT{j,i}(Press_Down_Rows{j,i},1);
        PRESS_UP_TIMES{j,i} = PRESS_STRUCT{j,i}(Press_Up_Rows{j,i},1);
        HOLD_Down_TIMES{j,i} = HOLD_STRUCT{j,i}(Hold_Down_Rows{j,i},1);
        HOLD_Up_TIMES{j,i} = HOLD_STRUCT{j,i}(Hold_Up_Rows{j,i},1);
        Press_num(j,i) = length(Press_Down_Rows{j,i});
        Hold_num(j,i) = length(Hold_Down_Rows{j,i});
        Press_Lever_Out_Time(j,i) = PRESS_STRUCT{j,i}(Forced_Press_Lever_Out_Rows(j,i),1);
        Hold_Lever_Out_Time(j,i) = HOLD_STRUCT{j,i}(Forced_Hold_Lever_Out_Rows(j,i),1);
    end
end

for i = 1:n
    for j = 1:5
        Hold_num1(j,i) = length(Hold_Up_Rows{j,i});
        Press_num(j,i) = length(Press_Up_Rows{j,i});
    end
end

for i = 1:n
    for j = 1:5
       if length(find(PRESS_STRUCT{j,i}(:,2) == 3330) == 1)
            PRESS_DURATIONS{j,i} = NaN;
       else    
            PRESS_DURATIONS{j,i} = PRESS_UP_TIMES{j,i} - PRESS_DOWN_TIMES{j,i};
       end
    end
end

for i = 1:n
    PRESS_DURATION_DATA{1,i} = [PRESS_DURATIONS{1,i};PRESS_DURATIONS{2,i};PRESS_DURATIONS{3,i};...
        PRESS_DURATIONS{4,i};PRESS_DURATIONS{5,i}];
end


for i = 1:n
   Mean_Press_Response_Durations(1,i) = nanmean(PRESS_DURATION_DATA{1,i}); 
end


for i = 1:n
   x_dur_pell{1,i} = 1:1:length(PRESS_DURATION_DATA{1,i}); 
end

for i = 1:n
    for j = 1:5
       if length(find(HOLD_STRUCT{j,i}(:,2) == 3330) == 1)
            HOLD_DURATIONS{j,i} = NaN;
       else    
            HOLD_DURATIONS{j,i} = HOLD_Up_TIMES{j,i} - HOLD_Down_TIMES{j,i};
       end
    end
end

for i = 1:n
    HOLD_DURATION_DATA{1,i} = [HOLD_DURATIONS{1,i};HOLD_DURATIONS{2,i};HOLD_DURATIONS{3,i};...
        HOLD_DURATIONS{4,i};HOLD_DURATIONS{5,i}];
end


for i = 1:n
   Mean_Hold_Response_Durations(1,i) = nanmean(HOLD_DURATION_DATA{1,i}); 
end


% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%- - Need Actual Press Up/Down Times to get IRT's
% - Three steps:
    % A) Getting the actual Times
    % B) Removing the first Down Time and Last UP Time
    % C) Difference between UP(n) and Down(n - 1)

    
% = = = For Presss = = =     
for i = 1:n
    for j = 1:5
       if length(find(PRESS_STRUCT{j,i}(:,2) == 3330) == 1)
            PRESS_FOR_IRT_UP_TIMES{j,i} = NaN;
            PRESS_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            PRESS_FOR_IRT_UP_TIMES{j,i} = PRESS_UP_TIMES{j,i};
            PRESS_FOR_IRT_DOWN_TIMES{j,i} = PRESS_DOWN_TIMES{j,i};
       end
    end
end

for i = 1:n
    for j = 1:5
       if length(find(PRESS_STRUCT{j,i}(:,2) == 3330) == 1)
            PRESS_FOR_IRT_UP_TIMES{j,i} = NaN;
            PRESS_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            PRESS_FOR_IRT_UP_TIMES{j,i}(end,:) = [];
            PRESS_FOR_IRT_DOWN_TIMES{j,i}(1,:) = [];
       end
    end
end


for i = 1:n
    for j = 1:5
       if length(find(PRESS_STRUCT{j,i}(:,2) == 3330) == 1)
            PRESS_IRTs{j,i} = NaN;
       else    
            PRESS_IRTs{j,i} = PRESS_FOR_IRT_DOWN_TIMES{j,i} - PRESS_FOR_IRT_UP_TIMES{j,i};
       end
    end
end

for i = 1:n
    PRESS_IRT_DATA{1,i} = [PRESS_IRTs{1,i};PRESS_IRTs{2,i};PRESS_IRTs{3,i};...
        PRESS_IRTs{4,i};PRESS_IRTs{5,i}];
end

for i = 1:n
   x_irt_pell{1,i} = 1:1:length(PRESS_IRT_DATA{1,i}); 
end

aa = 1:1:9;
aa = aa';
x_irt_num = [aa;aa;aa;aa;aa];



% = = = For Hold = = =     
for i = 1:n
    for j = 1:5
       if length(find(HOLD_STRUCT{j,i}(:,2) == 3330) == 1)
            HOLD_FOR_IRT_UP_TIMES{j,i} = NaN;
            HOLD_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            HOLD_FOR_IRT_UP_TIMES{j,i} = HOLD_Up_TIMES{j,i};
            HOLD_FOR_IRT_DOWN_TIMES{j,i} = HOLD_Down_TIMES{j,i};
       end
    end
end

for i = 1:n
    for j = 1:5
       if length(find(HOLD_STRUCT{j,i}(:,2) == 3330) == 1)
            HOLD_FOR_IRT_UP_TIMES{j,i} = NaN;
            HOLD_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            HOLD_FOR_IRT_UP_TIMES{j,i}(end,:) = [];
            HOLD_FOR_IRT_DOWN_TIMES{j,i}(1,:) = [];
       end
    end
end


for i = 1:n
    for j = 1:5
       if length(find(HOLD_STRUCT{j,i}(:,2) == 3330) == 1)
            HOLD_IRTs{j,i} = NaN;
       else    
            HOLD_IRTs{j,i} = HOLD_FOR_IRT_DOWN_TIMES{j,i} - HOLD_FOR_IRT_UP_TIMES{j,i};
       end
    end
end

for i = 1:n
    HOLD_IRT_DATA{1,i} = [HOLD_IRTs{1,i};HOLD_IRTs{2,i};HOLD_IRTs{3,i};...
        HOLD_IRTs{4,i};HOLD_IRTs{5,i}];
end

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% - - Finding Pauses and Bout Lengths


% = = = = NUMBER OF PAUSES AND PAUSE DURATIONS = = = = = = = 

% - Structure of Pause Rows for Different Durations
% - Row Storage Schedule: 1sec, 2sec, 5sec, 10sec
for i = 1:n
    Press_Pause_rows{1,i} = find(PRESS_IRT_DATA{1,i}(:,1) >= 1);
    Hold_Pause_rows{1,i} = find(HOLD_IRT_DATA{1,i}(:,1) >= 1);
    
    Press_Pause_rows{2,i} = find(PRESS_IRT_DATA{1,i}(:,1) >= 2);
    Hold_Pause_rows{2,i} = find(HOLD_IRT_DATA{1,i}(:,1) >= 2);
    
    Press_Pause_rows{3,i} = find(PRESS_IRT_DATA{1,i}(:,1) >= 5);
    Hold_Pause_rows{3,i} = find(HOLD_IRT_DATA{1,i}(:,1) >= 5);
    
    Press_Pause_rows{4,i} = find(PRESS_IRT_DATA{1,i}(:,1) >= 10);
    Hold_Pause_rows{4,i} = find(HOLD_IRT_DATA{1,i}(:,1) >= 10);
end

% =========================================================================================
% =========================================================================================
% = = Getting IRT's after removing the pauses 

Filtered_Press_IRT_DATA = {PRESS_IRT_DATA;PRESS_IRT_DATA;...
    PRESS_IRT_DATA;PRESS_IRT_DATA};
Filtered_Hold_IRT_DATA = {HOLD_IRT_DATA;HOLD_IRT_DATA;...
    HOLD_IRT_DATA;HOLD_IRT_DATA};

% = Getting IRT Data and Actually Removing the Different Pause Def's
for i = 1:n
    for j = 1:4
        if isempty(Press_Pause_rows{j,i})
            Filtered_Press_IRT_DATA{j,1}{1,i} = Filtered_Press_IRT_DATA{j,1}{1,i};
        else
            Filtered_Press_IRT_DATA{j,1}{1,i}(Press_Pause_rows{j,i}(:,1),:) = [];
        end
    end
end

for i = 1:n
    for j = 1:4
        if isempty(Hold_Pause_rows{j,i})
            Filtered_Hold_IRT_DATA{j,1}{1,i} = Filtered_Hold_IRT_DATA{j,1}{1,i};
        else
            Filtered_Hold_IRT_DATA{j,1}{1,i}(Hold_Pause_rows{j,i}(:,1),:) = [];
        end
    end
end


for i = 1:n
    for j = 1:4
        Mean_Press_IRT(j,i) = nanmean(Filtered_Press_IRT_DATA{j,1}{1,i});
    end
end

for i = 1:n
    for j = 1:4
        Mean_Hold_IRT(j,i) = nanmean(Filtered_Hold_IRT_DATA{j,1}{1,i});
    end
end
% =========================================================================================
% =========================================================================================



% - Determining the number of Pauses made
% - Row Storage Schedule: 1sec, 2sec, 5sec, 10sec
for i = 1:n
    for j = 1:4
        Press_Pause_num{j,i} = length(Press_Pause_rows{j,i});
        Hold_Pause_num{j,i} = length(Hold_Pause_rows{j,i});
    end
end

% - - Creating a Variable of just the "Pause" Durations
for i = 1:n
    for j = 1:4
        Press_Pause_Durations{j,i} = PRESS_IRT_DATA{1,i}(Press_Pause_rows{j,i},1);
        Hold_Pause_Durations{j,i} = HOLD_IRT_DATA{1,i}(Hold_Pause_rows{j,i},1);
    end
end

for i = 1:n
    for j = 1:4
        Mean_Press_Pause_Duration{j,i} = nanmean(Press_Pause_Durations{j,i});
        Mean_Hold_Pause_Duration{j,i} = nanmean(Hold_Pause_Durations{j,i});
        Sum_Press_Pause_Duration{j,i} = sum(Press_Pause_Durations{j,i});
        Sum_Hold_Pause_Duration{j,i} = sum(Hold_Pause_Durations{j,i});
    end
end

% = = = = = = Bout Lengths = = = = = = = = = 
% - Row Storage Schedule: 1sec, 2sec, 5sec, 10sec
for i = 1:n
    for j = 1:5

            Trial_Press_Pause_rows{j,i}{1,1} = find(PRESS_IRTs{j,i}(:,1) >= 1);
            Trial_Hold_Pause_rows{j,i}{1,1} = find(HOLD_IRTs{j,i}(:,1) >= 1);
            
            Trial_Press_Pause_rows{j,i}{2,1} = find(PRESS_IRTs{j,i}(:,1) >= 2);
            Trial_Hold_Pause_rows{j,i}{2,1} = find(HOLD_IRTs{j,i}(:,1) >= 2);
            
            Trial_Press_Pause_rows{j,i}{3,1} = find(PRESS_IRTs{j,i}(:,1) >= 5);
            Trial_Hold_Pause_rows{j,i}{3,1} = find(HOLD_IRTs{j,i}(:,1) >= 5);
            
            Trial_Press_Pause_rows{j,i}{4,1} = find(PRESS_IRTs{j,i}(:,1) >= 10);
            Trial_Hold_Pause_rows{j,i}{4,1} = find(HOLD_IRTs{j,i}(:,1) >= 10);
            
            Pel_No_Pause_Detector{j,i} = find(PRESS_IRTs{j,i}(:,1) >= .00001);
    end
end


% = = = TWO STEPS
% - A) Getting Opt Out Trials to Contain NaN's
% - B) Getting Trials with no Pauses to be able to lead to proper
% calculation of the bout lengths

% - A)
% - Getting trials which were optouts to contain NaN's
for i = 1:n
    for j = 1:5
        for k = 1:4
            if isnan(PRESS_DURATIONS{j,i})
                Trial_Press_Pause_rows{j,i}{k,1} = NaN;
            else
                Trial_Press_Pause_rows{j,i}{k,1} = Trial_Press_Pause_rows{j,i}{k,1};
            end
        end
    end
end

%- B) Getting Trials with no pauses to be able to calculate correct bout
%lengths

% - Number of Pauses per Trial
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Press_Pause_rows{j,i}{k,1})
                trial_num_press_pause{j,i}{k,1} = NaN;
            else
                trial_num_press_pause{j,i}{k,1} = length(Trial_Press_Pause_rows{j,i}{k,1});
            end
        end
    end
end


for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(HOLD_DURATIONS{j,i})
                trial_num_hold_pause{j,i}{k,1} = NaN;
            else
                trial_num_hold_pause{j,i}{k,1} = length(Trial_Hold_Pause_rows{j,i}{k,1});
            end
        end
    end
end

for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Hold_Pause_rows{j,i}{k,1})
                trial_num_hold_pause{j,i}{k,1} = NaN;
            else
                trial_num_hold_pause{j,i}{k,1} = length(Trial_Hold_Pause_rows{j,i}{k,1});
            end
        end
    end
end

% - - Preparing the data so I can extract the bout lengths from it
% - - 2 Step Process
        % A) Adding 1 to the pause rows so I can just subtract them and get the
               % bout lengths
        % B) Adding the total number of responses to the end 
        
% - - Getting the number needed for part B
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(PRESS_DURATIONS{j,i})
                Num_Press_Presses{j,i} = NaN;
            else
                Num_Press_Presses{j,i} = length(PRESS_DURATIONS{j,i});
            end
         end
     end
end
 
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(HOLD_DURATIONS{j,i})
                Num_Hold_Presses{j,i} = NaN;
            else
                Num_Hold_Presses{j,i} = length(HOLD_DURATIONS{j,i});
            end
         end
     end
 end

% - Implementing Part A and B
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Num_Press_Presses{j,i})
                Press_bout_length_trials_PRE{j,i}{k,1} = NaN;
            else
                Press_bout_length_trials_PRE{j,i}{k,1} = [1;Trial_Press_Pause_rows{j,i}{k,1};Num_Press_Presses{j,i}];
            end
        end
    end
end

for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Hold_Pause_rows{j,i}{k,1})
                Hold_bout_length_trials_PRE{j,i}{k,1} = NaN;
            else
                Hold_bout_length_trials_PRE{j,i}{k,1} = [1;Trial_Hold_Pause_rows{j,i}{k,1};Num_Hold_Presses{j,i}];
            end
        end
    end
end


% - - - Determining the length of the bouts within each trial
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Press_Pause_rows{j,i}{k,1})
                Press_Bout_Length_Trials{j,i}{k,1} = NaN;
            else
                Press_Bout_Length_Trials{j,i}{k,1} = diff(Press_bout_length_trials_PRE{j,i}{k,1});
            end
        end
    end
end

for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Hold_Pause_rows{j,i}{k,1})
                Hold_Bout_Length_Trials{j,i}{k,1} = NaN;
            else
                Hold_Bout_Length_Trials{j,i}{k,1} = diff(Hold_bout_length_trials_PRE{j,i}{k,1});
            end
        end
    end
end

% - Making variable with all pellet bout legnths over all trials
for i = 1:n
    for k = 1:4
        PRESS_BOUT_LENGTH_DATA{k,i} = [Press_Bout_Length_Trials{1,i}{k,1};...
            Press_Bout_Length_Trials{2,i}{k,1};...
            Press_Bout_Length_Trials{3,i}{k,1};...
            Press_Bout_Length_Trials{4,i}{k,1};...
            Press_Bout_Length_Trials{5,i}{k,1}];
    end
end



% - - Calculating the average pellet bout length

% = = = 
% = = =   Should Say Press = pel and Hold = suc but I was lazy
for i = 1:n
        Var_BOUT_1sec{1,i} = PRESS_BOUT_LENGTH_DATA{1,i};
        Var_BOUT_2sec{1,i} = PRESS_BOUT_LENGTH_DATA{2,i};
        Var_BOUT_5sec{1,i} = PRESS_BOUT_LENGTH_DATA{3,i};
        Var_BOUT_10sec{1,i} = PRESS_BOUT_LENGTH_DATA{4,i};
end

[Pel_BOUT_DF_1sec] = trials_cellmat(Var_BOUT_1sec,Var);
[Pel_BOUT_DF_2sec] = trials_cellmat(Var_BOUT_2sec,Var);
[Pel_BOUT_DF_5sec] = trials_cellmat(Var_BOUT_5sec,Var);
[Pel_BOUT_DF_10sec] = trials_cellmat(Var_BOUT_10sec,Var);


pb_DF_1sec = Pel_BOUT_DF_1sec.Trials_Matrix;
pb_DF_2sec = Pel_BOUT_DF_2sec.Trials_Matrix;
pb_DF_5sec = Pel_BOUT_DF_5sec.Trials_Matrix;
pb_DF_10sec = Pel_BOUT_DF_10sec.Trials_Matrix;


for i = 1:n
    Pel_NUM_BOUT(1,i) = length(pb_DF_1sec(:,i));
    Pel_NUM_BOUT(2,i) = length(pb_DF_2sec(:,i));
    Pel_NUM_BOUT(3,i) = length(pb_DF_5sec(:,i));
    Pel_NUM_BOUT(4,i) = length(pb_DF_10sec(:,i));
end

for i = 1:n
    Bout_NaN{1,i} = length(find(isnan(pb_DF_1sec(:,i)) == 1));
    Bout_NaN{2,i} = length(find(isnan(pb_DF_2sec(:,i)) == 1));
    Bout_NaN{3,i} = length(find(isnan(pb_DF_5sec(:,i)) == 1));
    Bout_NaN{4,i} = length(find(isnan(pb_DF_10sec(:,i)) == 1));
end

for i = 1:n
    Press_Number_of_BOUTS(1,i) = Pel_NUM_BOUT(1,i) - Bout_NaN{1,i};
    Press_Number_of_BOUTS(2,i) = Pel_NUM_BOUT(2,i) - Bout_NaN{2,i};
    Press_Number_of_BOUTS(3,i) = Pel_NUM_BOUT(3,i) - Bout_NaN{3,i};
    Press_Number_of_BOUTS(4,i) = Pel_NUM_BOUT(4,i) - Bout_NaN{4,i};
end

Mean_Press_Bout_Length(1,:) = nanmean(pb_DF_1sec);
Mean_Press_Bout_Length(2,:) = nanmean(pb_DF_2sec);
Mean_Press_Bout_Length(3,:) = nanmean(pb_DF_5sec);
Mean_Press_Bout_Length(4,:) = nanmean(pb_DF_10sec);


% - Making a vriable with all sucrose bout lengths combined
% - Making a vriable with all sucrose bout lengths combined
for i = 1:n
    for k = 1:4
        HOLD_BOUT_LENGTH_DATA{k,i} = [Hold_Bout_Length_Trials{1,i}{k,1};...
            Hold_Bout_Length_Trials{2,i}{k,1};...
            Hold_Bout_Length_Trials{3,i}{k,1};...
            Hold_Bout_Length_Trials{4,i}{k,1};...
            Hold_Bout_Length_Trials{5,i}{k,1}];
    end
end

% - Determining the average sucrose bout length
for i = 1:n
    Var_BOUT_suc_1sec{1,i} = HOLD_BOUT_LENGTH_DATA{1,i};
    Var_BOUT_suc_2sec{1,i} = HOLD_BOUT_LENGTH_DATA{2,i};
    Var_BOUT_suc_5sec{1,i} = HOLD_BOUT_LENGTH_DATA{3,i};
    Var_BOUT_suc_10sec{1,i} = HOLD_BOUT_LENGTH_DATA{4,i};
end

[Suc_BOUT_DF_1sec] = trials_cellmat(Var_BOUT_suc_1sec,Var);
[Suc_BOUT_DF_2sec] = trials_cellmat(Var_BOUT_suc_2sec,Var);
[Suc_BOUT_DF_5sec] = trials_cellmat(Var_BOUT_suc_5sec,Var);
[Suc_BOUT_DF_10sec] = trials_cellmat(Var_BOUT_suc_10sec,Var);

sb_DF_1sec = Suc_BOUT_DF_1sec.Trials_Matrix;
sb_DF_2sec = Suc_BOUT_DF_2sec.Trials_Matrix;
sb_DF_5sec = Suc_BOUT_DF_5sec.Trials_Matrix;
sb_DF_10sec = Suc_BOUT_DF_10sec.Trials_Matrix;

for i = 1:n
    Suc_NUM_BOUT(1,i) = length(sb_DF_1sec(:,i));
    Suc_NUM_BOUT(2,i) = length(sb_DF_2sec(:,i));
    Suc_NUM_BOUT(3,i) = length(sb_DF_5sec(:,i));
    Suc_NUM_BOUT(4,i) = length(sb_DF_10sec(:,i));
end

for i = 1:n
    Bout_NaN_suc{1,i} = length(find(isnan(sb_DF_1sec(:,i)) == 1));
    Bout_NaN_suc{2,i} = length(find(isnan(sb_DF_2sec(:,i)) == 1));
    Bout_NaN_suc{3,i} = length(find(isnan(sb_DF_5sec(:,i)) == 1));
    Bout_NaN_suc{4,i} = length(find(isnan(sb_DF_10sec(:,i)) == 1));
end

for i = 1:n
    Hold_Number_of_BOUTS(1,i) = Suc_NUM_BOUT(1,i) - Bout_NaN_suc{1,i};
    Hold_Number_of_BOUTS(2,i) = Suc_NUM_BOUT(2,i) - Bout_NaN_suc{2,i};
    Hold_Number_of_BOUTS(3,i) = Suc_NUM_BOUT(3,i) - Bout_NaN_suc{3,i};
    Hold_Number_of_BOUTS(4,i) = Suc_NUM_BOUT(4,i) - Bout_NaN_suc{4,i};
end

Mean_Hold_Bout_Length(1,:) = nanmean(sb_DF_1sec);
Mean_Hold_Bout_Length(2,:) = nanmean(sb_DF_2sec);
Mean_Hold_Bout_Length(3,:) = nanmean(sb_DF_5sec);
Mean_Hold_Bout_Length(4,:) = nanmean(sb_DF_10sec);




% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
for i = 1:n
    for j = 1:5
        if length(find(PRESS_STRUCT{j,i}(:,2) == 22)) == 1
            PRESS_Reward_Rows(j,i) = find(PRESS_STRUCT{j,i}(:,2) == 22);
        
        else
            PRESS_Reward_Rows(j,i) = NaN;
        end
    end
end

for i = 1:n
    for j = 1:5
        if length(find(HOLD_STRUCT{j,i}(:,2) == 25)) == 1
            HOLD_Reward_Rows(j,i) = find(HOLD_STRUCT{j,i}(:,2) == 25);
        else
            HOLD_Reward_Rows(j,i) = NaN;
        end
    end
end


for i = 1:n
    for j = 1:5
        if length(find(PRESS_STRUCT{j,i}(:,2) == 22)) == 1
            PRESS_Reward_Times(j,i) = PRESS_STRUCT{j,i}(PRESS_Reward_Rows(j,i),1); 
        else
            PRESS_Reward_Times(j,i) = NaN;
        end
    end
end

for i = 1:n
    for j = 1:5
        if length(find(HOLD_STRUCT{j,i}(:,2) == 25)) == 1
            HOLD_Reward_Times(j,i) = HOLD_STRUCT{j,i}(HOLD_Reward_Rows(j,i),1);
        else
            HOLD_Reward_Times(j,i) = NaN;
        end
    end
end



% - - Determining Times to reward from:
% - - 1) First press to reward
% - - 2) Trial start to reward

for i = 1:n
   for j = 1:5
      if isempty(PRESS_DOWN_TIMES{j,i}) == 1 
          Working_Press_Time(j,i) =  NaN;
          Total_Press_Time(j,i) =  NaN;
      else
          Working_Press_Time(j,i) =  PRESS_Reward_Times(j,i) - PRESS_DOWN_TIMES{j,i}(1,1);
          Total_Press_Time(j,i) =  PRESS_Reward_Times(j,i) - Press_Lever_Out_Time(j,i);
      end
   end
end

for i = 1:n
   for j = 1:5
       if isempty(HOLD_Down_TIMES{j,i}) == 1
           Working_Hold_Time(j,i) =  NaN;
           Total_Hold_Time(j,i) =  NaN;
       else
            Working_Hold_Time(j,i) =  HOLD_Reward_Times(j,i) - HOLD_Down_TIMES{j,i}(1,1);
            Total_Hold_Time(j,i) =  HOLD_Reward_Times(j,i) - Hold_Lever_Out_Time(j,i);
       end
   end
end


% - - Determining Latency to start working

for i = 1:n
   for j = 1:5
        if isempty(PRESS_DOWN_TIMES{j,i}) == 1 
              Latency_to_First_Press(j,i) = NaN;
        else
              Latency_to_First_Press(j,i) = PRESS_DOWN_TIMES{j,i}(1,1) - Press_Lever_Out_Time(j,i);
        end
   end
end

for i = 1:n
   for j = 1:5
      if isempty(HOLD_Down_TIMES{j,i}) == 1
            Latency_to_First_Hold(j,i) =  NaN;
      else    
            Latency_to_First_Hold(j,i) =  HOLD_Down_TIMES{j,i}(1,1) - Hold_Lever_Out_Time(j,i);
      end
   end
end


for i = 1:n
   Sub_PRESS_Reward_Times{1,i} = PRESS_Reward_Times(:,i); 
   Sub_HOLD_Reward_Times{1,i} = HOLD_Reward_Times(:,i);
end

SPRT = Sub_PRESS_Reward_Times;
SHRT = Sub_HOLD_Reward_Times;

pn = 5;
a = 1:1:pn;
a = a';

for i = 1:n
    Full_Press_Trials{1,i} = find(Press_num(:,i) == 20);
end

fpt = Full_Press_Trials;

for i = 1:n
    lfp(1,i) = length(Full_Press_Trials{1,i});
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - Getting Averages over Data

for i = 1:16
    Mean_Working_Press_Time = nanmean(Working_Press_Time);
    Mean_Total_Press_Time = nanmean(Total_Press_Time);
    Mean_Latency_to_First_Press = nanmean(Latency_to_First_Press);
    
    Mean_Working_Hold_Time = nanmean(Working_Hold_Time);
    Mean_Total_Hold_Time = nanmean(Total_Hold_Time);
    Mean_Latency_to_First_Hold = nanmean(Latency_to_First_Hold);
end

% - - - - - - - - - - DATA OUTPUT OF THE FUNCTION - - - - - - - - - - - %


DATA = struct('PRESS_DURATIONS_DATA',{{}},...
    'HOLD_DURATIONS_DATA',{{}},...
    'Mean_Press_DURATIONS',{{}},...
    'Mean_Hold_DURATIONS',{{}},...
    'Press_IRTs',{{}},...
    'Hold_IRTs',{{}},...
    'Mean_Press_IRTs',{{}},...
    'Mean_Hold_IRTs',{{}},...
    'Press_Pause_Durations',{{}},...
    'Hold_Pause_Durations',{{}},...
    'Press_Pause_num',{{}},...
    'Hold_Pause_num',{{}},...
    'Mean_Press_Pause_Duration',{{}},...
    'Mean_Hold_Pause_Duration',{{}},...
    'Sum_Press_Pause_Duration',{{}},...
    'Sum_Hold_Pause_Duration',{{}},...
    'trial_num_press_pause',{{}},...
    'trial_num_hold_pause',{{}},...
    'Press_Bout_Length_Trials',{{}},...
    'Hold_Bout_Length_Trials',{{}},...
    'PRESS_BOUT_LENGTH_DATA',{{}},...
    'HOLD_BOUT_LENGTH_DATA',{{}},...
    'Mean_Press_Bout_Length',{{}},...
    'Mean_Hold_Bout_Length',{{}}); 


DATA.PRESS_DURATIONS_DATA = PRESS_DURATION_DATA;   %Press Response Duations
DATA.HOLD_DURATIONS_DATA	= HOLD_DURATION_DATA; %Hold Response Durations

DATA.Mean_Press_DURATIONS = Mean_Press_Response_Durations;   %Mean response duration
DATA.Mean_Hold_DURATIONS	= Mean_Hold_Response_Durations;

DATA.Press_IRTs = Filtered_Press_IRT_DATA; %List of All Press IRT's
DATA.Hold_IRTs = Filtered_Hold_IRT_DATA;
DATA.Mean_Press_IRTs = Mean_Press_IRT; %Mean of All IRT's
DATA.Mean_Hold_IRTs = Mean_Hold_IRT;

DATA.Press_Pause_Durations	= Press_Pause_Durations; %Trial x trial - Pause Durations
DATA.Hold_Pause_Durations = Hold_Pause_Durations;

DATA.Press_Pause_num = Press_Pause_num;
DATA.Hold_Pause_num = Hold_Pause_num;

DATA.Mean_Press_Pause_Duration = Mean_Press_Pause_Duration; %Trial x Trial
DATA.Mean_Hold_Pause_Duration = Mean_Hold_Pause_Duration;	%Trial x Trial
DATA.Sum_Press_Pause_Duration = Sum_Press_Pause_Duration;	%Trial x Trial
DATA.Sum_Hold_Pause_Duration = Sum_Hold_Pause_Duration;	%Trial x Trial


DATA.trial_num_press_pause	= trial_num_press_pause; %Trial x Trial # Pauses
DATA.trial_num_hold_pause = trial_num_hold_pause;	%Trial x Trial # Pauses

DATA.Press_Bout_Length_Trials = Press_Bout_Length_Trials;	%Trial x Trial bout length		
DATA.Hold_Bout_Length_Trials = Hold_Bout_Length_Trials;	%Trial x Trial bout length


DATA.PRESS_BOUT_LENGTH_DATA = PRESS_BOUT_LENGTH_DATA;	%All bout lengths combined
DATA.HOLD_BOUT_LENGTH_DATA = HOLD_BOUT_LENGTH_DATA;	%%All bout lengths combined


DATA.Mean_Press_Bout_Length =	Mean_Press_Bout_Length;	% Averaged over all trials
DATA.Mean_Hold_Bout_Length = Mean_Hold_Bout_Length;

end


