function Cap_Stall_deletept(idx)

global Data
answer = questdlg(['Do you want to delete number ' num2str(idx) ' Caplillary?']);
if strcmp(answer,'Yes')
%     Data.Cap(idx,:) = [];
    if isfield(Data,'pts')
        if size(Data.pts,1) >= idx
            Data.pts(idx,:) = [];
        end
    end
    
%     if isfield(Data,pts1)
%         if size(Data.pts1,1) >= idx
%             Data.pts1(idx,:) = [];
%         end
%     end
%     if isfield(Data,pts2)
%         if size(Data.pts2,1) >= idx
%             Data.pts2(idx,:) = [];
%         end
%     end
%     if isfield(Data,pts3)
%         if size(Data.pts3,1) >= idx
%             Data.pts3(idx,:) = [];
%         end
%     end
    % delete column for that capillary in the autoStallingMatrix
    if isfield(Data,'AutoStallingMatrix')
        if size(Data.AutoStallingMatrix,1) >= idx
            Data.AutoStallingMatrix(idx,:) = [];
        end
    end
    
    % delete column for that capillary in the StallingMatrix
    if isfield(Data,'StallingMatrix')
        if size(Data.StallingMatrix,1) >= idx
            Data.StallingMatrix(idx,:) = [];
        end
    end
    
    if isfield(Data,'seg')
        if size(Data.seg,2) >= idx
            Data.seg(:,idx) = [];
        end
    end
    
    if isfield(Data,'Cap')
        if size(Data.Cap,1) >= idx
            Data.Cap(idx,:) = [];
        end
    end
    
    if isfield(Data,'Int_ts')
        if size(Data.Int_ts,1) >= idx
            Data.Int_ts(idx,:) = [];
        end
    end
    
    if isfield(Data,'GTStallingMatrix')
        if size(Data.GTStallingMatrix,1) >= idx
            Data.GTStallingMatrix(idx,:) = [];
        end
    end
    
    if isfield(Data,'ValidationFlag')
        if size(Data.ValidationFlag,1) >= idx
            Data.ValidationFlag(idx,:) = [];
        end
    end
    
    if isfield(Data,'segAnalysis')
        Data.segAnalysis.AveIntensity(idx,:) = [];
        Data.segAnalysis.AveCOV(idx,:) = [];
        Data.segAnalysis.DatasetAveInt(idx,:) = [];
        Data.segAnalysis.DatasetAveCOV(idx,:) = [];
    end
    
    handles = Data.handles;
    hObject = Data.hObject;
    eventdata = Data.eventdata;
    capStall('draw',hObject,eventdata,handles);
    if isfield(Data,'sliderobject')
        uicontrol(Data.sliderobject);
    end
end