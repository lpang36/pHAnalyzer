function [KOHpH,CaOH2pH] = myfunction (myvar)
%read from excel
[num,txt,rawKOH] = xlsread('raw.xlsx','data','B11:D11'); % koh room	=B3:F7 koh cold	=H3:L7 caoh2 room=N3:R7 caoh2 cold	=T3:X7
[num,txt,rawCaOH2] = xlsread('raw.xlsx','data','B12:D12'); %pls enter filename and range manually
[num,txt,rawCtrlKOH] = xlsread('raw.xlsx','data','B11:D11'); %pls enter filename and range manually, no blank cells 
[num,txt,rawCtrlCaOH2] = xlsread('raw.xlsx','data','B12:D12'); %pls enter filename and range manually, no blank cells 
[num,txt,rawCtrlPH] = xlsread('raw.xlsx','data', 'H10:H23'); % ph control strip
%convert to lch
KOH = hex2lch(rawKOH)
CaOH2 = hex2lch(rawCaOH2)
CtrlKOH = hex2lch(rawCtrlKOH)
CtrlCaOH2 = hex2lch(rawCtrlCaOH2)
CtrlPH = hex2lch (rawCtrlPH)
%avg controls
KOHavg = mean(CtrlKOH,3);
CaOH2avg = mean(CtrlCaOH2,3);
%compute distance - results returned in same format as excel cells, NaNs indicate empty cell
KOHresult = deltae(KOH,CtrlPH);
CaOH2result = deltae(CaOH2,CtrlPH);
%compute pH
KOHpH = interpolate(KOHresult)
CaOH2pH = interpolate(CaOH2result)


%subfunction for conversion
function output = hex2lch(input)
output = zeros(size(input));
for i = 1:size(input,1)
	for j = 1:size(input,2)
		if length(input(i,j))==1
			s = input(i,j);
            if(isnumeric(cell2mat(s)))
                s=num2str(cell2mat(s));
            else
                s=char(s);
            end
            s(1:2);
            s(3:4);
            s(5:6);
			r = hex2dec(s(1:2));
			g = hex2dec(s(3:4));
			b = hex2dec(s(5:6));
			lab = rgb2lab([r,g,b]);
			output(i,j,1) = lab(1);
			output(i,j,2) = sqrt(lab(2).^2+lab(3).^2);
			output(i,j,3) = atan2(lab(3),lab(2));
		else
			output(i,j,1) = 0;
            output(i,j,2) = 0;
            output(i,j,3) = 0;

		end
	end
end
%subfunction for distance metric
function output = deltae(data,CtrlPH)
output = zeros(size(data,1),size(data,2),14);

for i = 1:size(data,1)
	for j = 1:size(data,2)
		if sum(data(i,j,:))~= 0
            data(i,j)
			for k = 1:14
                %output(i,j,k) = sqrt(power(CtrlPH(k,1,1)-data(i,j,1),2)...
                %+power(CtrlPH(k,1,2)-data(i,j,2),2)...
                %+power(CtrlPH(k,1,3)-data(i,j,3),2));
                %^^ euclidean distance
                
                %using 1:1 
                %Calculate Sl
                if (CtrlPH(k,1,1)<16)
                    Sl = 0.511;
                else
                    Sl = 0.040975*CtrlPH(k,1,1)/(1+0.01765*CtrlPH(k,1,1));
                end
                
                %Calculate Sc
                Sc = 0.638 + (0.0638*CtrlPH(k,1,2)/(1+0.0131*CtrlPH(k,1,2)));
                
                %Calculate F
                F = sqrt(power(CtrlPH(k,1,2),4)/(power(CtrlPH(k,1,2),4)+1900));
                
                %Calculate T
                if (CtrlPH(k,1,3)<= 345 && CtrlPH(k,1,3)>= 164)
                    T = 0.56 + abs(0.2*cos(CtrlPH(k,1,3)+168));
                else
                    T = 0.36 + abs(0.4*cos(CtrlPH(k,1,3)+35));
                end
                
                %Calculate Sh
                Sh = Sc* (F*T +1 - F);
                
                %Calculate a's
                a1 = CtrlPH(k,1,2) * cos(CtrlPH(k,1,3));
                a2 = data(i,j,2) * cos(data(i,j,3));
                Da = a1-a2;
                %Calculate b's
                b1 = CtrlPH(k,1,2) * sin(CtrlPH(k,1,3));
                b2 = CtrlPH(k,1,2) * sin(CtrlPH(k,1,3));
                Db = b1-b2;
                
                %Calculate DC
                DC = CtrlPH(k,1,2) - data(i,j,2);
                
                %Calculate Delta H
                DH = sqrt(power(Da,2) + power(Db,2) + power(DC,2));
                
                %Final Equation
                output(i,j,k) = sqrt(power((data(i,j,1)-CtrlPH(k,1,1))/Sl ,2))...
                + power(( data(i,j,2)-CtrlPH(k,1,2))/Sc,2) ...
                + power((DH)/Sh,2);
                
            end 
		else
			output(i,j) = NaN;
		end
	end
end
function output = interpolate(data)
output = zeros(size(data,1),size(data,2))

for i = 1:size(data,1)
	for j = 1:size(data,2)
        %find min
        temp = 1;
        for k = 2:14
            if (data(i,j,k)<data(i,j,temp))
                temp = k;
            end
        end
        Min = data(i,j,temp);
        %find closest and interpolate
        try 
        if(abs(data(i,j,temp+1)-data(i,j,temp))<abs(data(i,j,temp)-data(i,j,temp-1)))
            %output(i,j) = temp-1/(temp-1+temp);
            output(i,j)  = temp +(data(i,j,temp)/(data(i,j,temp+1)+data(i,j,temp)));
        end
        if(abs(data(i,j,temp+1)-data(i,j,temp))>abs(data(i,j,temp)-data(i,j,temp-1)))
            %output(i,j) = temp+1/(temp+1+temp);
            output(i,j)  = temp -1 +(data(i,j,temp-1)/(data(i,j,temp-1)+data(i,j,temp)));
        end
        catch ME
            if temp == 1
                output(i,j)  = temp +(data(i,j,temp)/(data(i,j,temp+1)+data(i,j,temp)));
            else if temp ==14
                output(i,j)  = temp -1 +(data(i,j,temp-1)/(data(i,j,temp-1)+data(i,j,temp)));
                else
                    output (i,j) = 0;
                end
            end
        end
        
    end
end

        
        
        
        
