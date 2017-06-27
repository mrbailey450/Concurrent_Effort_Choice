% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% - This first line is to select the folder I want
folder_name = uigetdir

% - Now that you've selected the folder, execute all these lines of code
A = folder_name;

pth = strcat(A,'\');
rawlist = getrawdata(pth,'randycode');
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

[c_rawlist] = lev_prob(rawlist.data);

Var = c_rawlist;
X = rawlist.program;
[EFFORT_DATA] = Con_Effort_Choice_Data1(Var,X);



[DAT] = Con_Effort_Choice_Data(Var);
TS = DAT.trial_structure;

[CEC_TEMPORAL_DATA] = CEC_Temportal_Dat(Var,TS);
[CEC_BOUT_PAUSE_DATA] = CEC_Bout_Pause_Dat(Var,TS);


