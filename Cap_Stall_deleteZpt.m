function Cap_Stall_deleteZpt(idx)

global Data
answer = questdlg(['Do you want to delete number ' num2str(idx) ' Caplillary?']);
if strcmp(answer,'Yes')
    Data.pts(idx,:) = [];
    if isfield(Data,'Cap')
        if size(Data.Cap,1) >= idx
            Data.Cap(idx,:) = [];
        end
    end
%     idx = floor((num-1)/3)+1;
%     if rem(num,3) == 1
%         Data.pts1(idx,:) = [0 0 0];
%     elseif rem(num,3) == 2
%         Data.pts2(idx,:) = [0 0 0];
%     elseif rem(num,3) == 0
%         Data.pts3(idx,:) = [0 0 0];
%     end
%     if isfield(Data,'deletedpts')
%         Data.deletedpts = [Data.deletedpts; num];
%     else
%         Data.deletedpts = num;
%     end
    handles = Data.handles;
    hObject = Data.hObject;
    eventdata = Data.eventdata;
    Cap_Stall('draw',hObject,eventdata,handles);
    if isfield(Data,'sliderobjectZ')
        uicontrol(Data.sliderobjectZ);
    end
end
