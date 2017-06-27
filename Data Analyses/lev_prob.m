function [c_rawlist] = lev_prob(Var1)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% ==============PROBLEMATIC LEVER PRESSES=================%

n = length(Var1);

%Number of Lever Presses%
for i = 1:n
    Prob_Hold(1,i) =  length(find(Var1{1,i}(:,2) == 1015));
    Prob_Hold(2,i) =  length(find(Var1{1,i}(:,2) == 1017)); 
    Prob_Hold(3,i) = length(find(Var1{1,i}(:,2) == 1016));
    Prob_Hold(4,i) = length(find(Var1{1,i}(:,2) == 1018));
end;

%Identifying problematic hold cases%
for i= 1:n;
    prob(1,i) = Prob_Hold(1,i) - Prob_Hold(2,i);
    prob(2,i) = Prob_Hold(3,i) - Prob_Hold(4,i);
end;

% *************** Left Lever ***************** %

Left_Prob = find(prob(1,:)>=1)
num_L_Probs = length(Left_Prob);
nlp = num_L_Probs;

if length(nlp)>= 1
    for i = 1:nlp
        e = Var1{1,Left_Prob(1,i)};
            down = find (e(:,2) == 1015); up = find(e(:,2) == 1017);
            e(down(find(diff(down)==1)+1),:) = [];
        Var1{1,Left_Prob(1,i)} = e;
    end
else print('No Left Problems')
end    

%Left Lever%
if length(nlp)>= 1
    for i = 1:nlp
        e = Var1{1,Left_Prob(1,i)};
            down = find (e(:,2) == 1015); up = find(e(:,2) == 1017);
            if length(down)~=length(up)
            e(length(e(:,1))+1,:) = [max(e(:,1)) 1017];
            end
         Var1{1,Left_Prob(1,i)} = e;
    end
end

      
 % *************** Right Lever ***************** %   
    
Right_Prob = find(prob(2,:)>=1)
num_R_Probs = length(Right_Prob);
nrp = num_R_Probs;  
    
    
 if length(nrp)>= 1
    for i = 1:nrp
        e = Var1{1,Right_Prob(1,i)};
            down = find (e(:,2) == 1016); up = find(e(:,2) == 1018);
            e(down(find(diff(down)==1)+1),:) = [];
        Var1{1,Right_Prob(1,i)} = e;
    end
else print('No Left Problems')
end    

%Right Lever%
if length(nrp)>= 1
    for i = 1:nrp
        e = Var1{1,Right_Prob(1,i)};
            down = find (e(:,2) == 1016); up = find(e(:,2) == 1018);
            if length(down)~=length(up)
            e(length(e(:,1))+1,:) = [max(e(:,1)) 1018];
            end
         Var1{1,Right_Prob(1,i)} = e;
    end
end


for i = 1:n
    Prob_Hold(1,i) =  length(find(Var1{1,i}(:,2) == 1015));
    Prob_Hold(2,i) =  length(find(Var1{1,i}(:,2) == 1017)); 
    Prob_Hold(3,i) = length(find(Var1{1,i}(:,2) == 1016));
    Prob_Hold(4,i) = length(find(Var1{1,i}(:,2) == 1018));
end;
   

    %Identifying problematic hold cases%
for i= 1:n;
        prob(1,i) = Prob_Hold(1,i) - Prob_Hold(2,i);
        prob(2,i) = Prob_Hold(3,i) - Prob_Hold(4,i);
end;

    
    c_rawlist = Var1;
end

