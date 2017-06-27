function [CEC] = Con_Effort_Choice_Data1(Var,X)
%UNTITLED3 Summary of this function goes here

% = = = = = = = = Inputs = = = = = = = = = = = = = = = = 
% - Var = c_rawlist - need this bc hold problems and HT's
% - X = rawlist.program - Need this for the program names
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

n = length(Var);

% - - - - - - - Total Session Summary Variables - - - - - - - - - - -  


% Number of Reinforcers%
for i = 1:n
    Total_Session_Reins(1,i) = length(find(Var{1,i}(:,2) == 25));
end;

%Total Session Duration%
for i = 1:n;
    start(1,i) = find(Var{1,i}(:,2) == 113);
    stop(1,i) = find(Var{1,i}(:,2) == 114);
    Total_Session_Duration(1,i) = (Var{1,i}(stop(1,i),1) - Var{1,i}(start(1,i),1))/60;
end;

% - Determining the Lever ID code - (i.e Left = Hold or Press)
for i = 1:n
   if isempty(strfind(X{1,i}{1,1}, 'LH')) == 1
       Lever_Code_Name{1,i} = 'Left_Press';
       Lever_Code_ID(1,i) = 1;
   else
       Lever_Code_Name{1,i} = 'Right_Press';
       Lever_Code_ID(1,i) = 2;
   end
end

% - Determining total number of Left and Right Presses in entire session
for i = 1:n
   Total_Session_Presses(1,i) =  length(find(Var{1,i}(:,2) == 1015));
   Total_Session_Presses(2,i) = length(find(Var{1,i}(:,2) == 1016));
end

% - Using Hold vs Press code to differentiate Presses and Holds
for i = 1:n
   if Lever_Code_ID(1,i) == 1
       Total_Session_Lever_Presses(1,i) = Total_Session_Presses(1,i);
       Total_Session_Lever_Holds(1,i) = Total_Session_Presses(2,i);
   else
       Total_Session_Lever_Presses(1,i) = Total_Session_Presses(2,i);
       Total_Session_Lever_Holds(1,i) = Total_Session_Presses(1,i);
   end
end

%
%
% Add in Session Missed Dips
%
%
% - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - -






% - - - - - - - Forced Trials Variables - - - - - - - - - - - - 

%Number of OpOuts
for i = 1:n
    OpOuts(1,i) = length(find(Var{1,i}(:,2) == 3330));
end

%Forced Reinforcers

Forced_Reinforcers = 10 - OpOuts;

% - Finding all choice rows in session
for i  = 1:n
    Choice_Rows{1,i} = find(Var{1,i}(:,2) == 3301);
end

% - Defining row in whole session when choices begin
for i = 1:n
    Session_Choice_Start(1,i) = Choice_Rows{1,i}(1,1);
end

% - Easier name to work with
scs = Session_Choice_Start;

% - Determining total Presses and Holds in Forced Trials
for i = 1:n
    if Lever_Code_ID(1,i) == 1
        Total_Forced_Presses(1,i) = length(find(Var{1,i}(1:scs(1,i),2) ==  1015));
        Total_Forced_Holds(1,i) =  length(find(Var{1,i}(1:scs(1,i),2) ==  1016));
    else
        Total_Forced_Presses(1,i) = length(find(Var{1,i}(1:scs(1,i),2) ==  1016));
        Total_Forced_Holds(1,i) =  length(find(Var{1,i}(1:scs(1,i),2) ==  1015));
    end
end


% - Duration to complete Forced Trials
    

    % - Hack where End of Forced is Start of Choice
for i = 1:n;
    Forced_Start(1,i) = find(Var{1,i}(:,2) == 113);
end;
    
    %-Easier Name
    FS = Forced_Start;

% - Getting all ITI start rows bc 10th one will be end of Forced
for i = 1:n
   ITI_Start_Rows{1,i} =  find(Var{1,i}(:,2) == 121);
end

% - Defining the end row of Forced by 10th ITI start
for i = 1:n
    Forced_End_Row(1,i) = ITI_Start_Rows{1,i}(10,1);
end

    %Easier Name
    FER = Forced_End_Row;

% - How Long the Forced Trials Lasted (in munutes)
for i = 1:n
    Forced_Trials_Duration(1,i) = (Var{1,i}(FER(1,i),1) - Var{1,i}(FS(1,i),1))/60;
end


% - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - -




% - - - - - - - - - - - Choice Trials Variables - - - - - - - - - - - - 

%Number of Choice Trials
for i = 1:n
    Total_Choice_Trials(1,i) = length(find(Var{1,i}(:,2) == 3301));
end;

% - Finding all choice rows in session
for i  = 1:n
    Choice_Rows{1,i} = find(Var{1,i}(:,2) == 3301);
end

% - Defining row in whole session when choices begin
for i = 1:n
    Session_Choice_Start(1,i) = Choice_Rows{1,i}(1,1);
end

% - Easier name to work with
scs = Session_Choice_Start;

% - Determining total Presses and Holds in Choice Trials
    % - Row where choice starts in the session
    Start_Choice = scs + 1;
    % - Easier name to work with
    sc = Start_Choice;

% - - - Determining the number of Holds Made in Choice trials
for i = 1:n
   Total_Holds_in_Choice_Trials(1,i) = length(find(Var{1,i}(sc(1,i):end,2) ==  100));
end

% - Proportion of Choice Trials Which Were Holds

Proportion_Hold_Choices = Total_Holds_in_Choice_Trials ./ Total_Choice_Trials;

% - Total Reinforcers in the Choice Trials
for i = 1:n
   Total_Choice_Reins(1,i) = length(find(Var{1,i}(sc(1,i):end,2) ==  25));
end

% - Total # of Lever Presses and Holds in the Choice trials
for i = 1:n
    if Lever_Code_ID(1,i) == 1
        Total_Choice_Presses(1,i) = length(find(Var{1,i}(sc(1,i):end,2) ==  1015));
        Total_Choice_Holds(1,i) =  length(find(Var{1,i}(sc(1,i):end,2) ==  1016));
    else
        Total_Choice_Presses(1,i) = length(find(Var{1,i}(sc(1,i):end,2) ==  1016));
        Total_Choice_Holds(1,i) =  length(find(Var{1,i}(sc(1,i):end,2) ==  1015));
    end
end


% - Duration to complete Choice Trials
for i = 1:n
   Choice_Trials_Duration(1,i) =  (Var{1,i}(stop(1,i),1) - Var{1,i}(sc(1,i),1))/60;
end

% NUmber of Failed Hold Attempts

Number_Failed_Holds = Total_Choice_Holds - Total_Holds_in_Choice_Trials;

% - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - -




% = = = = = = = = = = = Outputs of the Function = = = = = = = = = = = = = 
CEC = struct('Session_Totals',{{}},'Forced_Trials',{{}},'Choice_Trials',{{}});

CEC.Session_Totals = struct('Total_Session_Reins',[],...
    'Total_Session_Duration',[],...
    'Total_Session_Lever_Presses',[],...
    'Total_Session_Lever_Holds',[]);
CEC.Forced_Trials = struct('OpOuts',[],...
    'Forced_Reinforcers',[],...
    'Total_Forced_Presses',[],...
    'Total_Forced_Holds',[],...
    'Forced_Trials_Duration',[]);
CEC.Choice_Trials = struct('Total_Holds_in_Choice_Trials',[],...
    'Total_Choice_Trials',[],...
    'Proportion_Hold_Choices',[],...
    'Total_Choice_Reins',[],...
    'Total_Choice_Presses',[],...
    'Total_Choice_Holds',[],...
    'Choice_Trials_Duration',[],...
    'Number_Failed_Holds',[]);

CEC.Session_Totals.Total_Session_Reins = Total_Session_Reins;
CEC.Session_Totals.Total_Session_Duration = Total_Session_Duration;
CEC.Session_Totals.Total_Session_Lever_Presses = Total_Session_Lever_Presses;
CEC.Session_Totals.Total_Session_Lever_Holds = Total_Session_Lever_Holds;

CEC.Forced_Trials.OpOuts = OpOuts; 
CEC.Forced_Trials.Forced_Reinforcers = Forced_Reinforcers;
CEC.Forced_Trials.Total_Forced_Presses = Total_Forced_Presses;
CEC.Forced_Trials.Total_Forced_Holds = Total_Forced_Holds;
CEC.Forced_Trials.Forced_Trials_Duration = Forced_Trials_Duration;

CEC.Choice_Trials.Total_Holds_in_Choice_Trials = Total_Holds_in_Choice_Trials;
CEC.Choice_Trials.Total_Choice_Trials = Total_Choice_Trials;
CEC.Choice_Trials.Proportion_Hold_Choices = Proportion_Hold_Choices;
CEC.Choice_Trials.Total_Choice_Reins = Total_Choice_Reins;
CEC.Choice_Trials.Total_Choice_Presses = Total_Choice_Presses;
CEC.Choice_Trials.Total_Choice_Holds = Total_Choice_Holds;
CEC.Choice_Trials.Choice_Trials_Duration = Choice_Trials_Duration;
CEC.Choice_Trials.Number_Failed_Holds = Number_Failed_Holds;
end

