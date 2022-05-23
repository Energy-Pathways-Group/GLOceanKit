function RepeatingAndNameValue(mynum, names,options)
arguments
    mynum (1,1) = 0
end
arguments (Repeating)
    names char
end
arguments
    options.allcaps (1,1) = 0
end

fprintf('mynum: %f\n',mynum);

for iVar=1:length(names)
    fprintf('name: %s\n',names{iVar})
end
if options.allcaps == 1
    fprintf('Should print in all caps\n');
else
    fprintf('keep case\n');
end

end