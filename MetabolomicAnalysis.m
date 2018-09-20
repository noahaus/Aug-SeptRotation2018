%author: Noah Legall
%This code is very specific so it can not be ran on a machine.
%my hope was to show that I do indeed have proficency with MATlab


%% Start of analysis (Control + Experimental Data Loading)
% Control Dataset. E07A was a failed inoculation and used as an example in this pipeline 
%(It wouldn't make sense to use this in real life) so the baseline for the
% targeted metabolites should reflect a normal state.
sExp = 'E07A';
sSubjects = 'All'; sSpecies = 'Mm';sTissue = 'Whole Blood';

%sProtocol = 'QM_FIA'; 
sProtocol = 'QM_LCMS';
disp('-----------------------------------------------')
disp('Get metabolomic data for all subjects')
strucData_cont = funGetProtocolData(cStructures,sExp,sSubjects,sProtocol,'sTissue',sTissue);

% I want to start with replacing all null values from the dataset
% A is that matrix
A = strucData_cont.cRawCounts{1,1};
%dimensions of the dataset
[m,n] = size(A);
%straightforward nested for loops to find nulls and replace them with zeros
for row = 1:m
    for column = 1:n
        if isnan(A(row,column))
           A(row,column) = 0;
        end
    end
end

% Experimental Dataset. E03 infected monkeys with malarial parasite and
% noticed changes in targeted metabolomics.
sExp = 'E03';
sSubjects = {'RZe13','RUn13','RTi13','RWr13'}; sSpecies = 'Mm';sTissue = 'Whole Blood';

%sProtocol = 'QM_FIA'; 
sProtocol = 'QM_LCMS';
disp('-----------------------------------------------')
disp('Get metabolomic data for all subjects')
strucData_exp = funGetProtocolData(cStructures,sExp,sSubjects,sProtocol,'sTissue',sTissue);

% I want to start with replacing null values from the dataset
% B is that matrix
B = strucData_exp.cRawCounts{1,1};
%dimensions of the dataset
[m,n] = size(B);
%straightforward nested for loops to find nulls and replace them with zeros
for row = 1:m
    for column = 1:n
        if isnan(B(row,column))
           B(row,column) = 0;
        end
    end
end

%% Analyze the Same Elements
% Lets start with comparing what are the targeted metabolites that are in
% both lists.


A_cell = num2cell(A);
control_cell = num2cell(strucData_cont.VarNames{1,1});
cont_cell = horzcat(control_cell,A_cell);

B_cell = num2cell(B);
experimental_cell = num2cell(strucData_exp.VarNames{1,1});
exp_cell = horzcat(experimental_cell,B_cell);

% which cell array is larger?

%if control is larger, then this is the dataset with the extra metabolites
if size(control_cell,1) > size(experimental_cell,1)
    %process the cell array to extract the string value of the first column
    tempCell = {};
    for i=1:size(exp_cell,1)
        tempCell =[tempCell; exp_cell{i,1}{1,1}]; 
    end
    
    
    [m,n] = size(cont_cell);
    %if the string element is in the smaller cell array, we keep the value
    %otherwise, delete the row that string was on
    for row = 1:m
        if ismember(cell2mat(cont_cell{row,1}),tempCell)
            continue
        else
            cont_cell(row,:) = [];
        end
    end
elseif size(control_cell,1) == size(experimental_cell,1)
%Same rationale if experimental dataset is larger than the control.
else
   tempCell = {};
    for i=1:size(cont_cell,1)
        tempCell =[tempCell; cont_cell{i,1}{1,1}]; 
    end 
    
    [m,n] = size(exp_cell);
    for row = 1:m
        if ismember(cell2mat(exp_cell{row,1}),tempCell)
            continue
        else
            exp_cell{row,:} = {};
        end
    end
end
%% Two Sample T test
% This is the statistical test to see which metabolites are enriched
% formula is this:
% t = (x1-x2)/sqrt([(s1^2/n1)+(s2^2/n2)])
% fortunately, there is already a t-test function in matlab

%convert the cell arrays into proper matricies 
control_raw_score = cell2mat(cont_cell(:,2:7));
experimental_raw_score = cell2mat(exp_cell(:,2:5));

%We want to not include the zeros to the analysis, so let's filter them out
%each time with an array
control_values = [];
experiment_values = [];

%The enriched array will be binary with 1 signifying enriched and 0 meaning
%not enriched
enriched = [];

[m,n] = size(control_raw_score);
[j,k] = size(experimental_raw_score);


for row = 1:m
    for column = 1:n
       if control_raw_score(row,column) == 0
           continue;
       else
           control_values = [control_values;control_raw_score(row,column)];
       end
    end
    
    for column = 1:k
       if experimental_raw_score(row,column) == 0
           continue;
       else
           experiment_values = [experiment_values;experimental_raw_score(row,column)];
       end
    end
    
    %This is where the ttest occurs. luckily MATlab has this function
    %already created
    [h,p] = ttest2(control_values,experiment_values,'Alpha',0.05);
    enriched = [enriched;h];
end

experimental_raw_score = horzcat(experimental_raw_score,enriched);
