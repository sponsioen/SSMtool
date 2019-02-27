
function hwait=updatewaitbar(hwait,i,n)
if nargout
    if ~ishandle(hwait)
        hwait = waitbar(0,'Don''t close me >:(');
    end
end
if ~mod(i,max(1,floor(n/100)))
    waitbar(i/n,hwait);
end