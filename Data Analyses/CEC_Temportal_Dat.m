function [CEC_TEMPORAL_DATA] = CEC_Temportal_Dat(Var,TS)
%UNTITLED2 Summary of this function goes here

%   Var = rawlist.data1
%   TS = CON_EFFORT_DATA.trial_structure
%   


%   Detailed explanation goes here

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

n = length(Var);

%LLeverOn = 0027 \***
%RLeverOn= 0028  \***

%HoldTrial = 3320
%PressTrial = 3321

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

% - - Determining the Press Rows and Hold Rows

for i = 1:n
    for j = 1:5
        if length(find(PRESS_STRUCT{j,i}(:,2) == 28)) == 0
            Press_Rows{j,i} = find(PRESS_STRUCT{j,i}(:,2) == 1015);
            Forced_Press_Lever_Out_Rows(j,i) = find(PRESS_STRUCT{j,i}(:,2) == 27);
        else
            Press_Rows{j,i} = find(PRESS_STRUCT{j,i}(:,2) == 1016);
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

% - Determining the number of lever presses and holds
for i = 1:n
    for j = 1:5
        Press_num(j,i) = length(Press_Rows{j,i});
        Hold_num(j,i) = length(Hold_Down_Rows{j,i});
    end
end

% - Getting the times of all the lever presses
for i = 1:n
    for j = 1:5
        if Press_num(j,i) == 0
            PRESS_TIMES{j,i} = NaN;
        else
            PRESS_TIMES{j,i} = PRESS_STRUCT{j,i}(Press_Rows{j,i},1);
        end
    end
end

% - Getting the times of all the lever holds
for i = 1:n
    for j = 1:5
        if Hold_num(j,i) == 0
            HOLD_Down_TIMES{j,i} = NaN;
            HOLD_Up_TIMES{j,i} = NaN;
        else
            HOLD_Down_TIMES{j,i} = HOLD_STRUCT{j,i}(Hold_Down_Rows{j,i},1);
            HOLD_Up_TIMES{j,i} = HOLD_STRUCT{j,i}(Hold_Up_Rows{j,i},1);
        end
    end
end


% - Getting the time the lever is inserted 
for i = 1:n
    for j = 1:5
        Press_Lever_Out_Time(j,i) = PRESS_STRUCT{j,i}(Forced_Press_Lever_Out_Rows(j,i),1);
        Hold_Lever_Out_Time(j,i) = HOLD_STRUCT{j,i}(Forced_Hold_Lever_Out_Rows(j,i),1);
    end
end

%HOLD_Down_TIMES{5,12}(end,:) = [];
% - Calculating the duration of each hold
for i = 1:n
    for j = 1:5
       HOLD_DURATIONS{j,i} = HOLD_Up_TIMES{j,i} - HOLD_Down_TIMES{j,i};
    end
end

for i = 1:n
    for j = 1:5
        if length(find(PRESS_STRUCT{j,i}(:,2) == 25)) == 1
            PRESS_Reward_Rows(j,i) = find(PRESS_STRUCT{j,i}(:,2) == 25);
        
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
    
    
% - - Finding the time that Press Trials get rewarded
for i = 1:n
    for j = 1:5
        if length(find(PRESS_STRUCT{j,i}(:,2) == 25)) == 1
            PRESS_Reward_Times(j,i) = PRESS_STRUCT{j,i}(PRESS_Reward_Rows(j,i),1); 
        else
            PRESS_Reward_Times(j,i) = NaN;
        end
    end
end

% - - Finding the time that Hold Trials get rewarded
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
      Working_Press_Time(j,i) =  PRESS_Reward_Times(j,i) - PRESS_TIMES{j,i}(1,1);
      Total_Press_Time(j,i) =  PRESS_Reward_Times(j,i) - Press_Lever_Out_Time(j,i);
      Working_Hold_Time(j,i) =  HOLD_Reward_Times(j,i) - HOLD_Down_TIMES{j,i}(1,1);
      Total_Hold_Time(j,i) =  HOLD_Reward_Times(j,i) - Hold_Lever_Out_Time(j,i);
   end
end


% - - Determining Latency to start working

for i = 1:n
   for j = 1:5
      Latency_to_First_Press(j,i) = PRESS_TIMES{j,i}(1,1) - Press_Lever_Out_Time(j,i);
      Latency_to_First_Hold(j,i) =  HOLD_Down_TIMES{j,i}(1,1) - Hold_Lever_Out_Time(j,i);
      Latency_to_Press_Reward(j,i) = PRESS_Reward_Times(j,i) - Press_Lever_Out_Time(j,i);
      Latency_to_Hold_Reward(j,i) = HOLD_Reward_Times(j,i) - Hold_Lever_Out_Time(j,i);
   end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - Getting Averages over Data

for i = 1:16
    Mean_Working_Press_Time = nanmean(Working_Press_Time);
    Mean_Total_Press_Time = nanmean(Total_Press_Time);
    Mean_Latency_to_First_Press = nanmean(Latency_to_First_Press);
    
    Mean_Working_Hold_Time = nanmean(Working_Hold_Time);
    Mean_Total_Hold_Time = nanmean(Total_Hold_Time);
    Mean_Latency_to_Hold_Press = nanmean(Latency_to_First_Hold);
end


% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% - - - - - - - - - - - Output of the Function - - - - - - - - - - - - 

CEC_TEMPORAL_DATA = struct('PRESS_TIMES',{{}},...
    'HOLD_Down_TIMES',{{}},...
    'HOLD_Up_TIMES',{{}},...
    'Press_num',[],...
    'Hold_num',[],...
    'HOLD_DURATIONS',{{}},...
    'PRESS_Reward_Times',[],...
    'HOLD_Reward_Times',[],...
    'Working_Press_Time',[],...
    'Total_Press_Time',[],...
    'Working_Hold_Time',[],...
    'Total_Hold_Time',[],...
    'Latency_to_First_Press',[],...
    'Latency_to_First_Hold',[],...
    'Data_Structs',{{}},...
    'Means',{{}});


CEC_TEMPORAL_DATA.PRESS_TIMES = PRESS_TIMES;
CEC_TEMPORAL_DATA.HOLD_Down_TIMES = HOLD_Down_TIMES;
CEC_TEMPORAL_DATA.HOLD_Up_TIMES = HOLD_Up_TIMES;
CEC_TEMPORAL_DATA.Press_num = Press_num;
CEC_TEMPORAL_DATA.Hold_num = Hold_num;

CEC_TEMPORAL_DATA.HOLD_DURATIONS = HOLD_DURATIONS;


CEC_TEMPORAL_DATA.PRESS_Reward_Times = PRESS_Reward_Times;
CEC_TEMPORAL_DATA.HOLD_Reward_Times = HOLD_Reward_Times;


CEC_TEMPORAL_DATA.Working_Press_Time = Working_Press_Time;
CEC_TEMPORAL_DATA.Total_Press_Time = Total_Press_Time;
CEC_TEMPORAL_DATA.Working_Hold_Time = Working_Hold_Time;
CEC_TEMPORAL_DATA.Total_Hold_Time = Total_Hold_Time;


CEC_TEMPORAL_DATA.Latency_to_First_Press = Latency_to_First_Press;
CEC_TEMPORAL_DATA.Latency_to_First_Hold = Latency_to_First_Hold;

% - - The Press and Hold Data Structs for Forced Trials - - - - 
CEC_TEMPORAL_DATA.Data_Structs = struct('PRESS_STRUCT',{{}},'HOLD_STRUCT',{{}});
CEC_TEMPORAL_DATA.Data_Structs.PRESS_STRUCT = PRESS_STRUCT;
CEC_TEMPORAL_DATA.Data_Structs.HOLD_STRUCT = HOLD_STRUCT;


CEC_TEMPORAL_DATA.Means = struct('Mean_Working_Press_Time',[],...
    'Mean_Total_Press_Time',[],...
    'Mean_Latency_to_First_Press',[],...
    'Mean_Working_Hold_Time',[],...
    'Mean_Total_Hold_Time',[],...
    'Mean_Latency_to_Hold_Press',[]);




CEC_TEMPORAL_DATA.Means.Mean_Working_Press_Time = Mean_Working_Press_Time;
CEC_TEMPORAL_DATA.Means.Mean_Total_Press_Time = Mean_Total_Press_Time;
CEC_TEMPORAL_DATA.Means.Mean_Latency_to_First_Press = Mean_Latency_to_First_Press;
    
CEC_TEMPORAL_DATA.Means.Mean_Working_Hold_Time = Mean_Working_Hold_Time;
CEC_TEMPORAL_DATA.Means.Mean_Total_Hold_Time = Mean_Total_Hold_Time;
CEC_TEMPORAL_DATA.Means.Mean_Latency_to_Hold_Press = Mean_Latency_to_Hold_Press;







end

