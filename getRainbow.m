function [ color ] = getRainbow( Npath )
% This function allows you to get an RGB vector for visualizing graphs, 
% depending on the minimum required number of colors - Npath.
color = [];
if Npath<=12
    step = 0.5;
elseif Npath<=18
    step = 1/3;
elseif Npath<=24
    step = 0.25;
elseif Npath<=30
    step = 0.2;
else
    step = 0.1;
end;

for r=0:step:(1-step)
    color = [color; [r,1,0]];
end;
for g=flip((0+step):step:1)
    color = [color; [1,g,0]];
end;
for b=0:step:(1-step)
    color = [color; [1,0,b]];
end;
for r=flip((0+step):step:1)
    color = [color; [r,0,1]];
end;
for g=0:step:(1-step)
    color = [color; [0,g,1]];
end;
for b=flip((0+step):step:1)
    color = [color; [0,1,b]];
end;
end

