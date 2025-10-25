%gcf = figure;
k=[];
set(gcf,'keypress','k=get(gcf,''currentchar'');');
while 1
  %do stuff here
  if ~isempty(k)
    if strcmp(k,'s'); break; end;
    if strcmp(k,'p'); pause; k=[]; end;
  end
end