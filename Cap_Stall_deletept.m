function Cap_Stall_deletept(idx)

global Data
answer = questdlg(['Do you want to delete number ' num2str(idx) ' Caplillary?']);
if strcmp(answer,'Yes')
    Data.Cap(idx,:) = [];
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
    handles = Data.handles;
    hObject = Data.hObject;
    eventdata = Data.eventdata;
    Cap_Stall('draw',hObject,eventdata,handles);
    if isfield(Data,'sliderobject')
        uicontrol(Data.sliderobject);
    end
end