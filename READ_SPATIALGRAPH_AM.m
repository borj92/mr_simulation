%READ_SPATIALGRAPH_AM.m 
%This program reads Amira Graph files 
clear all
clc
%Opening the Amira file
%str1 = 'Amira_Graph_Example.txt';
%str1 = 'Flow2Amira(Small).txt';
str1 = input('Enter filename of amira data:');
fileID = fopen(str1);
curLine = fgetl(fileID);
mgc = '# AmiraMesh 3D ASCII 2.0';
if strcmp(mgc,curLine) == 0
    error_Message = 'Error, READ_SPAITALGRAPH_AM: Not a valid Amira file';
    disp(error_Message)
    fclose(fileID);
    return 
end 
fclose(fileID);

%Reading the header
fileID = fopen(str1);
curLine = fgetl(fileID);
i = 1;
j = 1;
vdefinition(j).vname = [];
 definition(i).type = [];
while ischar(curLine)
    
    %Defintion Line
    k = defCheck(curLine);
    if isempty(k)
    else
        spl = strsplit(curLine,' ');
        nspl = numel(spl);
        if nspl ~= 3
            disp('Error, READ_SPATIALGRAPH_AM: Unexpected formatting on defintion line');
            return
        else
            definition(i).type = spl(2);
            definition(i).len = spl(3);
            i = i+1;
        end
    end

    %Parameter Line
    if isempty(paramCheck(curLine))
    else
        paramLine = fgetl(fileID);
        pspl = strsplit(paramLine,' ');
        npspl = numel(pspl);
        mgc = '"HxSpatialGraph"';
        if strcmp(pspl(npspl),mgc) == 0
            disp('Error, READ_SPATIALGRAPH_AM: Not a SpatialGraph file!')
        end
    end
 
    %Variable Definition line
    if isempty(atCheck(curLine))
    else
        vspl = strsplit(curLine,' ');
        nvspl = numel(vspl);
        if nvspl == 6
            vdefinition(j).vtype = vspl{1};
            tmp = strsplit(vspl{3},{'[',']'});
            ntmp = numel(tmp);
            if ntmp == 3
                vdefinition(j).vdatatype = tmp{1};
                vdefinition(j).vdatadim = str2num(tmp{2});
            else
                vdefinition(j).vdatatype = tmp{1};
                vdefinition(j).vdatadim = 1;
            end
            vdefinition(j).vname = vspl{4};
            vdefinition(j).vmarker = vspl{6};
        end
        j= j+1;
    end
    
    curLine = strtrim(fgetl(fileID));
    %End of header line
    if isempty(hashCheck(curLine))
    else 
        curLine = 0; 
    end 

end
fclose(fileID);

for l = 1:j-1
    if isempty(vdefinition(l).vname)
    disp('Error, READ_SPATIALGRAPH_AM: No variables defined in header')
    end 
end  

%Reading the data 
fileID = fopen(str1);
line = fgetl(fileID);
dvals2 = [];
dvals3 = [];
dvals4 = [];
dvals5 = [];
dvals6 = [];
dvals7 = [];
while ischar(line)
    if atCheck(line) == 1
        i = 1;
        checkspl = strsplit(line,'@');
        n = str2num(checkspl{2});
        data(n).Name = vdefinition(n).vname;
        data(n).Marker = vdefinition(n).vmarker;
        data(n).Dim = vdefinition(n).vdatadim;
        dline = strtrim(fgetl(fileID));
        while isempty(atCheck(dline)) && ischar(dline)
            if dline == -1
                break
            end
            dval = strsplit(dline,' ');
            if isempty(dval)
            else
                switch n
                    case 1
                        if isempty(unicode2native(dval{1}))
                        else
                            for j = 1:numel(dval);
                                dvals(i,j) = str2num(dval{j});
                            end
                        end
             
                    case 2
                        if isempty(unicode2native(dval{1}))
                        else 
                            for j = 1:2
                                dvals2(i,j) = str2num(dval{j});
                            end
                        end
                        
                    case 3
                        a = numel(dval);
                        bytes = unicode2native(dval{1});
                        if isempty(bytes)
                        else
                            dvals3(i,1) = str2num(dval{1});
                        end
                        
                    case 4
                        bytes = unicode2native(dval{1});
                        if isempty(bytes)
                        else
                            for j = 1:numel(dval)
                                dvals4(i,j) = str2num(dval{j});
                            end
                        end
                        
                    case 5
                        bytes = unicode2native(dval{1});
                        if isempty(bytes)
                        else
                            dvals5(i,1) = str2num(dval{1});
                        end
                        
                    case 6
                        bytes = unicode2native(dval{1});
                        if isempty(bytes)
                        else
                            for j = 1:numel(dval)
                                dvals6(i,j) = str2num(dval{j});
                            end
                        end
                    case 7
                        bytes = unicode2native(dval{1});
                        if isempty(bytes)
                        else
                            for j = 1:numel(dval)
                                dvals7(i,j) = str2num(dval{j});
                            end
                        end
                end 
                data(1).Val = dvals;
                data(1).NumEl = numel(dvals);
                data(2).Val = dvals2;
                data(2).NumEl = numel(dvals2);
                data(3).Val = dvals3;
                data(3).NumEl = numel(dvals3);
                data(4).Val = dvals4;
                data(4).NumEl = numel(dvals4);
                data(5).Val = dvals5;
                data(5).NumEl = numel(dvals5);
                if isempty(dvals6)
                else
                    data(6).Val = dvals6;
                    data(6).NumEl = numel(dvals6);
                end
                if isempty(dvals7)
                else
                    data(7).Val = dvals7;
                    data(7).NumEl = numel(dvals7);
                end

                dline = fgetl(fileID);
                if dline == -1
                    break
                end
                dline = strtrim(dline);
                i = i+1;
            end
        end
        status = fseek(fileID, -7,0);
        disp('Done with stage')
        disp(n)
    end
    line = fgetl(fileID);
end
for i = 1:5
    data(i).NumEl =  data(i).NumEl/data(i).Dim;
end
fclose(fileID);
%save Amira_Data.mat data definition vdefinition -v7.3
save Amira_DataS2.mat data definition vdefinition -v7.3
clear all
load ('Amira_DataS2.mat')
disp('File Read: Complete')