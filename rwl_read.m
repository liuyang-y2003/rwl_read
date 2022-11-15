function [TRm,yr,cores,log]=rwl_read(varargin)
% Read decadal tree-ring measurement data (.rwl) and correct formating errors 
% <Syntax>
% [TRm,yr,cores,log]=rwl_read(filename)
% [TRm,yr,cores,log]=rwl_read(___,'round')
% [TRm,yr,cores,log]=rwl_read(___,'zero')
% 
% <Input>
% filename - Name and extension of the file to import, specified as a 
%     character vector or a string scalar.
%
% <Optional Input>
% 'round' - Option to round decimals, specified as 'filled'.
% 'zero' - Option to treat zeros as missing values, specified as 'zero'.
% 
% <Output>
% TRm - Tree-ring measurement data from the file, returned as an n-by-p matrix, 
%     where each row of TRm represents one year, and each column represents
%     one core.
% yr - Years, returned as a vector whose size is n.
% cores - Core ids, returned as a cell array whose size is p.
% log - Error log, returned as a two-column matrix. The first colume is the
%     error type and the second column is the line number where the error occurred.
%     error type list
%     1. Time error
%     2. Core ID error
%     3. No data row
%     4. Duplicate core
%     5. Stop marker absence
%     6. Special characters
%     7. Unexpected missing value
%     8. Segmented core series
%
% <Dependencies>
% This function requires insertarray.m.
% 
% <Examples>
% [TRm,yr,cores,log]=rwl_read('test.rwl');
% [TRm,yr,cores,log]=rwl_read('test.rwl','round','zero');
% 
% <Copyright>
% Author:   Yang Liu
% Contact:  liu_yang@igsnrr.ac.cn
% Update:   2022-11-15
% Version:  1.0.1

log = []; % set default
filename = varargin{1};
[rd,z0] = deal(false); % round decimals (rd), treat 0 as missing value (z0)

k = 2;
while k<=length(varargin)
    switch lower(varargin{k})
        case 'round'
            rd = true;
        case 'zero'
            z0 = true;
    end
    k = k+1;
end

fileID = fopen(filename); % Read .rwl file
Dcell = textscan(fileID,'%s','delimiter',{'\n','\r\n','\r'},'MultipleDelimsAsOne',1,'whitespace',''); % Data in cell format
Dcell = Dcell{:};
fclose(fileID);

line = (1:length(Dcell))'; % line numbers used for log
Dcell = deblank(Dcell);

% errpr 3.1: No data row, remove lines with no ring data (length<=12),
% e.g. test.rwl, Line 85; mong007t.rwl, Line 113
ind = cellfun(@(x) length(x), Dcell)<=12;
if sum(ind)>0
    log = [log;[ones(sum(ind),1)*3.1,line(ind)]];
    Dcell(ind)=[];
    line(ind)=[];
end

Dchar = char(Dcell); % Data in char format

% cheak header lines
if all(Dchar(1:3,8) == ['1';'2';'3']) % Regular 3 header lines
    Dchar(1:3,:)=[];
    line(1:3)=[];
end

% errpr 3.2: No data row, remove notation lines in mid rows, including unregularly header lines
% e.g. test.rwl, Line 34; swed303.rwl, Line 77
ind = sum(isletter(Dchar(:,9:72)),2)>0;
if sum(ind)>0
    temp = [ones(sum(ind),1)*3.2,line(ind)];
    log = [log;temp(temp(:,2)>3,:)];
    Dchar(ind,:)=[];
    line(ind)=[];
end

% deal with years earlier than -1000 (occupy 5 columns), e.g. ca714.rwl, Line 354
if ismember('-',Dchar(:,8))
    ind = Dchar(:,8)=='-';
    temp = Dchar(:,8);
    temp(ind)=' ';
    Ccell= cellstr([Dchar(:,1:7),temp]); 
    Dchar(~ind,8)=' '; % in case of both name and minus sign in col 8
    yr1=8; % year starts from col 8
else
    Ccell= cellstr(Dchar(:,1:8));
    yr1=9; % year starts from col 9
end

% error 6.1: Special characters, replace special chars (non digit, '.', '-', and space) with space
% e.g. test.rwl, Line 78; ital026.rwl, Line 117
temp = Dchar(:,yr1:72);
ind = ~isstrprop(temp,'digit') & temp~='.' & temp~=' ' & temp~='-';
if sum(ind(:))>0
    temp(ind) = ' ';
    Dchar(:,yr1:72)=temp;
    ind = sum(ind,2)>0;
    log = [log;[ones(sum(ind),1)*6.1,line(ind)]];
end

% error 6.2: Special characters, replace '.' with out adjacent digits by space
% e.g. test.rwl, Line 79; al002l.rwl, Line 325
temp = Dchar(:,13:end);
ind = int8((temp=='.'));
ind(:,1:end-1) = ind(:,1:end-1) + int8(temp(:,2:end)==' '); % check the column in the right
ind(:,end) = ind(:,end)*2; % in case of . in the end
ind = ind>=2 & temp=='.';
if sum(ind(:))>0
    temp(ind) = ' ';
    Dchar(:,13:end)=temp;
    ind = sum(ind,2)>0;
    log = [log;[ones(sum(ind),1)*6.2,line(ind)]];
end

% check line # for a new core by id names
% 'case insensitive', e.g.: ph001e.rwl, Lines 58, 59, 76 and 93
new_id = [true;~strcmpi(Ccell(1:end-1),Ccell(2:end))];
% check line # for a new core by years, in which the decades is not continuous with the previous line
new_yr = [true;floor(str2num(Dchar(1:end-1,yr1:12))/10) ~= floor(str2num(Dchar(2:end,yr1:12))/10)-1];
% check line # for a new core by stop marker of 999 or -9999
[~,new_99] = regexp(cellstr(Dchar(:,13:72)),'( 999|-9999)\s*$'); % end with ' 999' or '-9999'
new_99 = cellfun(@isempty,new_99);
new_99 = [true;~new_99(1:end-1)];

% error 2: Core ID error
% 2.1; detect and correct core ID shifts
% e.g. test.rwl, Lines 9-10 and 17-18; brit095.rwl, Line 132
ind = find(new_id & ~new_yr & ~new_99); % only id change
new_yr99 = new_yr & new_99;
while ~isempty(ind)
    Cb = find(new_id(1:ind(1)-1),1,'last'); % same core back
    Cf = find(new_id(ind(1):end),2)+ind(1)-1; % same core formard
    if length(Cf) == 2 % in case ind(1) is the last line
        Cf = Cf(2)-1;
    else
        Cf = Cf(1);
    end
    dvdL = find(new_yr99(Cb+1:Cf))+Cb-1; % actual divided line
    if length(dvdL) == 1
        Ccell(Cb:dvdL)=Ccell(Cb);
        Ccell(dvdL+1:Cf)=Ccell(Cf);
        if ind(1)>dvdL
            range = unique([ind(1)-1,dvdL+1]);
        else
            range = unique([ind(1),dvdL]);
        end
        log = [log;[ones(range(end)-range(1)+1,1)*2.1,line(range(1):range(end))]];
        new_id = [true;~strcmpi(Ccell(1:end-1),Ccell(2:end))];
        ind = find(new_id & ~new_yr & ~new_99);
    else
        ind(1)=[];
    end
end
% 2.2: detect and correct simple core ID error
% e.g. test.rwl, Line 25; ca654.rwl, Line 751
ind = find(new_id & ~new_yr & ~new_99); % only id change
while ~isempty(ind)
    Ccell(ind(1)) = Ccell(ind(1)-1);
    log = [log;[2.2,line(ind(1))]];
    new_id = [true;~strcmpi(Ccell(1:end-1),Ccell(2:end))];
    ind = find(new_id & ~new_yr & ~new_99);
end

% error 1: Time error
ind = find([new_yr(1:end-1) & new_yr(2:end);new_yr(end)]);
new_id99 = [new_id & new_99;true];
new3 = new_id | new_yr | new_99;
if ~isempty(ind)
    for i = 1:length(ind)
        % current core's middle year is wrong, e.g. test.rwl, Line 42; cana015.rwl, Line 533
        if ind(i)>1 && ind(i)< size(Dchar,1) && ~new_id99(ind(i)) && ~new_id99(ind(i)+1) ...
                && floor(str2num(Dchar(ind(i)-1,yr1:12))/10)+2 == floor(str2num(Dchar(ind(i)+1,yr1:12))/10)
            new_yr(ind(i):ind(i)+1)=false;
            yrr = num2str(str2num(Dchar(ind(i)+1,yr1:12))-10);
            Dchar(ind(i),13-length(yrr):12)=yrr;
            log = [log;[1.1,line(ind(i))]];
        end
        % a new core's first year is wrong, e.g. test.rwl, Line 35; mexi043e.rwl, Line 925
        if ind(i)< size(Dchar,1)-1 && new_id99(ind(i)) && ~new_id99(ind(i)+1) && ~new3(ind(i)+2) % check if it has only one row, 
            yrr = num2str(str2num(Dchar(ind(i)+1,yr1:12))-10);
            rem = Dchar(ind(i),12); % maintain ones place
            Dchar(ind(i),13-length(yrr):12)=yrr;
            Dchar(ind(i),12)=rem;
            new_yr(ind(i)+1)=0;
            log = [log;[1.2,line(ind(i))]];
        end
        % current core's last year is wrong, e.g. test.rwl, Line 49; ak104.rwl, Lines 1328 and 1365
        if ind(i)>2 && ~new_id99(ind(i)) && new_id99(ind(i)+1) && ~new3(ind(i)-1)
            yrr = num2str(str2num(Dchar(ind(i)-1,yr1:12))+10);
            Dchar(ind(i),13-length(yrr):12)=yrr;
            Dchar(ind(i),12)='0'; % in case of the previous row is the first row
            new_yr(ind(i))=0;
            log = [log;[1.3,line(ind(i))]];
        end
    end
end
new3 = new_id | new_yr | new_99; % rows that a new core stats

cores = Ccell(new3); %  col cell of continuous unique names
new3 = cumsum(new3);

Ncores = length(cores); % number of separate cores

% Compute the first and last row of each core, and the number of rows
Crows = [find([1;diff(new3)]),diff([find([1;diff(new3)]);length(new3)+1])]; 

yr = unique(floor(str2num(Dchar(:,yr1:12))/10)*10);
yr = (yr(1):yr(end)+9)';
TRm = nan(length(yr),Ncores);

% check length of consecutive digit columns, it there is no space between two values, insert a space column
% e.g. eth003.rwl, line 9
temp = movsum(Dchar(:,13:72)~=' ',6,2);
% this needs Image Toolbox, do not use
% temp = bwconncomp(insertarray(1,Dchar(:,13:72),repmat(' ',size(Dchar,1)-1,60),1.5:size(Dchar,1))~=' ');
% temp = cellfun(@(x) length(x),temp.PixelIdxList);
if max(temp(:))>5
    Dchar = insertarray(2,Dchar,repmat(' ',size(Dchar,1),10),12.5:6:66.5);
    L1 = 7;
else
    L1 = 6;
end

% read data
for n=1:Ncores % loop over cores
    nrows = Crows(n,2) ;  % number of decade lines for this core
    rowstart = Crows(n,1); % first row of this series in C
    X = [];
    for m = 1:nrows % loop over decades
        if m == nrows || m == 1
            c = deblank(Dchar(rowstart+m-1,1:12+L1*10)); % deblank the first/last row
        else
            c = Dchar(rowstart+m-1,1:12+L1*10); % don't deblank middle rows,e.g. russ1.rwl, Line 728, a missing value (blank) in the end
        end
        rem = mod(str2double(c(yr1:12)),10); % get ones place
        
        % some files use consecutive spaces for missing values
        % e.g. test.rwl, Line 53; cana590.rwl, Line 5329
        gap = strfind(c(13:end),repmat(' ',1,L1)); % L1-digit spaces
        if ~isempty(gap) && length(str2num(c(13:end))) < length(c(13:end))/L1 % have blank
            log = [log;[7.1,line(rowstart+m-1)]];
            if rem > 0 && m ~=1 % if the row is not the first row but the year is not a full decade, bring blank to the left and data to the right
                % e.g. japa017.rwl, Line 749
                c = [c(1:12),c(13+(10-rem)*L1:end),c(13:12+(10-rem)*L1)];
                gap = strfind(c(13:end),repmat(' ',1,L1));
            end
            while ~isempty(gap)
                c(9+L1+gap(1):11+L1+gap(1)) = 'NaN'; % replace blank gaps by NaN
                gap = strfind(c(13:end),repmat(' ',1,L1));
            end
        end
        temp = str2num(c(13:end))';
        
        if m == 1
            % remove leading 9990, e.g. brit5.rwl, Line 4
            if temp(1) == 9990 && find(diff(temp),1) == rem
                temp(1:rem)=[];
                log = [log;[7.2,line(rowstart+m-1)]];
            end
            
            % if there is a gap in the end of the first row, insert NaN to match year
            % e.g. japa017.rwl, Line 748
            if (temp(end) == 999 || temp(end) == -9999)
                temp(end) = [];
                temp = [temp;nan(10-length(temp)-rem,1)];
            end
            
            % error 1.4: check if the start year matchs the number of data in the first row
            % e.g. test.rwl, Line 76; az560.rwl, Line 111
            if rem ~= 10 - length(temp)
                yr1st = str2double(c(yr1:12)) - rem +10 - length(temp);
                log = [log;[1.4,line(rowstart+m-1)]];
            else
                yr1st=str2double(c(yr1:12));
            end
        end
        
        if m == nrows
            % remove trailing 9990, e.g. brit5.rwl, Line 19
            if temp(end) == 9990
                log = [log;[7.2,line(rowstart+m-1)]];
                temp(find(diff(temp),1,'last')+1:end)=[];
            end
            
            if temp(end) == 999 || temp(end) == -9999
                temp(end) = [];
            else if m ~= 1 && isempty(strfind(Dchar(rowstart+m-1,12+L1*10:end),'999')) % sometimes the stop marker is after column 72, e.g. russ032w.rwl, line 412
                    log = [log;[5,line(rowstart+m-1)]];
                end
            end
            if m == 1
                yrend = str2double(c(yr1:12)) + length(temp) - 1;
            else
                yrend = str2double(c(yr1:12)) - mod(str2double(c(yr1:12)),10) + length(temp) - 1;
            end
        end
        X = [X;temp];
    end
    TRm((yr1st:yrend)-yr(1)+1,n) = X;
end
cores=deblank(cores);

ind = find(sum(~isnan(TRm),2),1);
if ind>1 % remove leading rows with all NaNs
    yr(1:ind-1)=[];
    TRm(1:ind-1,:)=[];
end
ind = find(sum(~isnan(TRm),2),1,'last');
if ind<length(yr) % remove trailing rows with all NaNs
    yr(ind+1:end)=[];
    TRm(ind+1:end,:)=[];
end

% error 8: Segmented core series
% e.g. test.rwl, Lines 60-64; indi002.rwl, Lines 200 and 212
seq = find([1;~strcmp(cores(1:end-1),cores(2:end))]);
seq = [seq,[seq(2:end)-1;length(cores)],diff([seq;length(cores)+1])];
seq = seq(seq(:,3)>1,:);
if ~isempty(seq)
    for i = size(seq,1):-1:1
        if all(sum(~isnan(TRm(:,seq(i,1):seq(i,2))),2)<2)
            log = [log;[ones(seq(i,3),1)*8,line(Crows(seq(i,1):seq(i,2),1))]];
            TRm(:,seq(i,1)) = mean(TRm(:,seq(i,1):seq(i,2)),2,'omitnan');
            TRm(:,seq(i,1)+1:seq(i,2))=[];
            cores(seq(i,1)+1:seq(i,2))=[];
            Crows(seq(i,1)+1:seq(i,2),:)=[];
        end
    end
end

% determine remaining 999 is real data value or missing value
if z0
    TRm(TRm<=0)=nan;
else
    TRm(TRm<0)=nan;
end
if ~isempty(find(TRm==999,1))
    if ~isempty(find((TRm>900 & TRm<1100 & TRm~=999),1))
        TRm(TRm==999) = nan; % change to NaN
    end
end

ind = sum(~isnan(TRm),1)==0; % remove cores with all NaN values
cores(ind)=[];
TRm(:,ind)=[];
Crows(ind,:)=[];
if isempty(TRm)
    log = [];
end

% error 4: Duplicate core
% e.g. test.rwl, Lines 24-28 and 55-59; mexi067.rwl, Line 350
[~,~,ic] = unique(cores);
TR2 = [ic';TRm]';
TR2(isnan(TR2)) = -0.1;
[~,ia,~] = unique(TR2,'stable','rows');
if length(ia)<size(TRm,2)
    log = [log;[ones(size(TRm,2)-length(ia),1)*4.1,line(Crows(setdiff(1:size(TRm,2),ia),1))]];
    cores = cores(ia);
    TRm = TRm(:,ia);
    Crows = Crows(ia,:);
end
% e.g. test.rwl, Lines 91-95
[~,ia,~] = unique(cores);
if length(ia)<size(TRm,2)
    log = [log;[ones(size(TRm,2)-length(ia),1)*4.2,line(Crows(setdiff(1:size(TRm,2),ia),1))]];
end

if ~isempty(log)
    log = floor(log);
end

if rd
    TRm = round(TRm);
end
