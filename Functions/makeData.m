function [Data,Im] = makeData(MacroTime,ArrTime,NLine,NColumn,TFrame,DT,Freq)
%makeData takes photon arrival times and generate the format used by BNP code.
%
%INPUT:
%   MacroTime: photon arrival times with respect to start of experiment (ms)
%   ArrTime: photon arrival times with respect to start of pulses (photonPhase) (ms)
%   NLine:   number of rows in the image data frame (scanning trajectories)
%   NColumn: number of columns in the image data frame
%   TFrame:  vector containing the scanning times of frames (frameTime) (ms)
%   Dt:      pixel duel times (time that a pixel was scanned) (ms)
%   Freq:    laser frequency
%
%OUTPUT:
%   Data:    Data in a format used by BNP, which is structure array where
%            every element correspond to a pixel. The fields are
%      Dt:     photon arrival times (with respect to start of pulse) 
%              associated to the pixel (ns).    
%      W:      Binary paramter indicating that a pulse was empty or not.
%   Im:      Image frame of the data where pixel values are number of
%            detected photons for pixels
%

Data(NLine,NColumn).Dt = [];
Data(NLine,NColumn).W = [];
NPulse = 0;
Im = zeros(1,NLine*NColumn);

for frame = 1:length(TFrame)-1
    %ID of photons from this frame
    ID = MacroTime<TFrame(frame+1) & MacroTime>TFrame(frame);
    %photons arrivals from this frame with respect to start of the experiment
    ThisFrameMicroTimes = MacroTime(ID) - TFrame(frame);
    %photon arrivals from this frame with respect to start of pulses
    ThisArrTime = ArrTime(ID); 
    %beginning and end of pixel scans 
    TBins = 0:DT:TFrame(frame+1)-TFrame(frame) + DT/2;
    %pixel ID of photons for this frame
    %[~,Ind] = histc(ThisFrameMicroTimes,TBins);
    [Values,~,Ind] = histcounts(ThisFrameMicroTimes,TBins);
    %H = histogram(ThisFrameMicroTimes,TBins);
    if length(Values) ~= NLine*NColumn
        %check if the scan was completed if not go to the next frame
        continue 
    end
    %Assigning photons to pixels
    tID = 0;
    for mm = 1:NColumn
        for nn = 1:NLine
            tID = tID + 1;
            Data(nn,mm).Dt = cat(1,Data(nn,mm).Dt,ThisArrTime(Ind==tID)'/Freq/255*10^9); 
            Data(nn,mm).X_Confocal = mm - 0.5;
            Data(nn,mm).Y_Confocal = nn - 0.5;
            Data(nn,mm).Z_Confocal = 0;
        end
    end
    Im = Im + Values;
    NPulse = NPulse + DT*10^(-3)*Freq;
end

for nn = 1:NLine*NColumn
    Data(nn).W = zeros(round(NPulse),1); 
    Data(nn).W(1:Im(nn)) = 1;
end
Im = reshape(Im,NLine,NColumn);
end